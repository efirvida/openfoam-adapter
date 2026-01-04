---
title: Configure the OpenFOAM adapter
permalink: adapter-openfoam-config.html
keywords: adapter, openfoam, configuration, preciceDict, controlDict
summary: "Write a system/preciceDict, set compatible boundary conditions, and activate the adapter in your system/controlDict."
---

In order to run a coupled simulation, you need to:

1. prepare a preCICE configuration file (described in the [preCICE configuration](https://precice.org/configuration-overview.html)),
2. prepare an adapter's configuration file,
3. set the coupling boundaries in the OpenFOAM case,
4. load the adapter, and
5. start all the solvers normally, from the same directory, e.g., in two different terminals.

You may skip the section _"Advanced configuration"_ in the beginning, as it only concerns special cases.

## The adapter's configuration file

The adapter is configured via the file `system/preciceDict`. This file is an OpenFOAM dictionary with the following form:

```c++
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      preciceDict;
}

preciceConfig "precice-config.xml";

participant Fluid;

modules (FSI);

interfaces
{
  Interface1
  {
    mesh              Fluid-Mesh;
    patches           (interface);
    locations         faceCenters;

    readData
    (
      Displacement
    );

    writeData
    (
      Force
    );
  };
};
```

The `participant` needs to be the same as the one specified in the `preciceConfig`,
which is the main preCICE configuration file. The `preciceConfig` can be a path and needs to be wrapped with quotation marks.

The list `modules` should contain `FSI` for fluid-structure interaction simulations.

In the `interfaces`, we specify the coupling interfaces (here only one).
The `mesh` needs to be the same as the one specified in the `preciceConfig`.
The `patches` specifies a list of the names of the OpenFOAM boundary patches that are
participating in the coupled simulation. These need to be defined in the files
included in the `0/` directory. The names of the interfaces (e.g., `Interface1`) are arbitrary and are not used.

The `locations` field is optional and its default value is `faceCenters` (with `faceCentres` also accepted), signifying that the interface mesh is defined on the cell face centers. An alternative option is `faceNodes`, which defines the mesh on the face nodes and is needed, e.g., for reading displacements in an FSI scenario.

For fluid-structure interaction, coupled quantities can be:

- `writeData`:
  - fluid participants: `Force`
  - solid participants: `Displacement`
- `readData`:
  - fluid participants: `Displacement`, `AngularVelocity`
  - solid participants: `Force`

## Configuration of the OpenFOAM case

A few changes are required in the configuration of an OpenFOAM case, in order to specify the interfaces and load the adapter. For some solvers, additional parameters may be needed (see "advanced configuration").

### Boundary conditions

The type of the `readData` needs to be compatible with the respective boundary
conditions set for each field in the `0/` directory of the case.

Read the [OpenFOAM User Guide](https://www.openfoam.com/documentation/user-guide/boundaries.php) for more on boundary conditions.

#### FSI

- For `readData(Displacement)`, you need the following:
  - `type movingWallVelocity` for the interface (e.g., `flap`) in `0/U`,
  - `type fixedValue` for the interface (e.g., `flap`) in the `0/pointDisplacement`, and
  - `solver displacementLaplacian` in the `constant/dynamicMeshDict`. The solver [`RBFMeshMotionSolver` from solids4foam is also known to work](https://github.com/precice/openfoam-adapter/pull/241), since the OpenFOAM adapter v1.2.0 and the solids4foam v2.0.

```c++
// File 0/U
interface
{
    type            movingWallVelocity;
    value           uniform (0 0 0);
}

// File 0/pointDisplacement
interface
{
    type            fixedValue;
    value           $internalField;
}

// File constant/dynamicMeshDict
dynamicFvMesh       dynamicMotionSolverFvMesh;
motionSolverLibs    ("libfvMotionSolvers.so");
solver              displacementLaplacian;
```

### Load the adapter

To load this adapter, you must include the following in
the `system/controlDict` configuration file of the case:

```c++
libs ("libpreciceAdapterFunctionObject.so");
functions
{
    preCICE_Adapter
    {
        type preciceAdapterFunctionObject;
        errors strict; // optional
    }
}
```

This directs the solver to use the `preciceAdapterFunctionObject` function object,
which is part of the `libpreciceAdapterFunctionObject.so` shared library.
The name `preCICE_Adapter` can be arbitrary.

The `errors strict` option is optional and [available since OpenFOAM v2012](https://www.openfoam.com/news/main-news/openfoam-v20-12/post-processing#post-processing-function-object-error-handling). Since the adapter is necessary to do a coupled simulation, this option instructs OpenFOAM to stop in case it faces issues with loading the adapter. For OpenFOAM versions that don't support this, remove the option.

If you are using other function objects in your simulation, add the preCICE adapter to the end of the list. The adapter will then be executed last, which is important, as the adapter also controls the end of the simulation. When the end of the simulation is detected, the adapter also triggers the `end()` method of all function objects.

***

## Advanced configuration

These additional parameters may only concern some users in special cases. Keep reading if you want to use the nearest-projection mapping, if you are using a solver with different variable names, or if you are trying to debug a simulation.

### Nearest-projection mapping

The [preCICE documentation](https://precice.org/couple-your-code-defining-mesh-connectivity.html) contains a detailed description of nearest-projection mappings in preCICE. In summary, we need to explicitly enable the `connectivity` option to create edges between the interface mesh points and give them to preCICE:

```c++
interfaces
{
  Interface1
  {
    mesh              Fluid-Mesh-Centers;
    locations         faceCenters;
    connectivity      false;
    patches           (interface);

    // ... writeData, readData ...
  };

  Interface2
  {
    mesh              Fluid-Mesh-Nodes;
    locations         faceNodes;
    connectivity      true;
    patches           (interface);

    // ... writeData, readData ...
  };
};
```

This `connectivity` boolean is optional and defaults to `false`. Note that `connectivity true` can only be used with `locations faceNodes`.

Even if the coupling data is associated to `faceCenters` in the solver, we can select `faceNodes` as locations type: the respective data will be interpolated from faces to nodes. Also, connectivity is only needed and supported for `writeData`. Therefore, we need to split the interface in a "read" and a "write" part, as shown above.

### Additional properties for FSI solvers

The adapter's FSI functionality supports both compressible and incompressible solvers, as well as solid (e.g., solids4Foam) solvers.

For incompressible solvers, it tries to read uniform values for the density and kinematic viscosity (if it is not already available) from the `FSI` subdictionary of `preciceDict`:

```c++
nu              nu [ 0 2 -1 0 0 0 0 ] 1e-03;
rho             rho [1 -3 0 0 0 0 0] 1;
```

#### User-defined solver type

The adapter tries to automatically determine the solver type,
based on the dictionaries that the solver uses.
However, you may manually specify the solver type to be `incompressible`, `compressible`, or `solid` (e.g., for solids4Foam) for an FSI simulation:

```c++
FSI
{
    solverType solid;
}
```

#### Parameters and fields with different names

The names of the parameters and fields that the adapter looks for
can be changed, in order to support a wider variety of solvers.
You may specify the following parameters in the adapter's configuration
file (the values correspond to the default values):

```c++
FSI
{
    // Displacement fields
    namePointDisplacement pointD;
    nameCellDisplacement D;
    // Force field
    nameForce Force; // For solids4Foam: solidForce
}
```

Use the option `namePointDisplacement unused;` for solvers that do not create a pointDisplacement field, such as the RBFMeshMotionSolver.

#### Restarting FSI simulations

Restarting a coupled simulation using the OpenFOAM adapter works in principle in the same way as restarting a standalone OpenFOAM solver. However, the adapter and preCICE define the coupling interface and the interface node locations based on the mesh at the particular time. In case of FSI simulations, the interface deforms over time, which leads to the definition of a deformed interface during the restart. The boolean variable `restartFromDeformed` (`true` by default) allows to account for the previously accumulated interface deformation such that the initial interface configuration (`t = 0`) is completely recovered. The setting here needs to agree with the behavior of the selected solid solver.

```c++
FSI
{
    // Account for previous displacement during a restart
    restartFromDeformed true;
}
```

{% important %}
The option here defines the way the interface mesh is initialized when restarting an FSI simulation in OpenFOAM. In order to restart a coupled simulation, your solid solver needs to be capable of restarting as well. Furthermore, the two participants need to follow the same assumption for the initialization, which for OpenFOAM you can configure with this option. You can find more information about restarting coupled simulations on [Discourse](https://precice.discourse.group/t/how-can-i-restart-a-coupled-simulation/675).
{% endimportant %}

#### Debugging

The user can toggle debug messages at [build time](https://precice.org/adapter-openfoam-get.html).

## Coupling OpenFOAM with 2D solvers

The adapter asks preCICE for the dimensions of the coupling data defined in the `precice-config.xml` (2D or 3D). It then automatically operates in either 3D (normal) or 2D (reduced) mode, with z-axis being the out-of-plane dimension. [Read more](https://github.com/precice/openfoam-adapter/pull/96). In 2D mode, the adapter also supports axisymmetric cases.

## Porting your older cases to the current configuration format

In earlier versions of the adapter, we were using a yaml-based configuration format,
with the adapter configuration file usually named as `precice-adapter-config.yml`.
We moved to a OpenFOAM dictionary format in [#105](https://github.com/precice/openfoam-adapter/pull/105),
to reduce the dependencies. You may also find the [tutorials #69](https://github.com/precice/tutorials/pull/69)
to be a useful reference (file changes).
