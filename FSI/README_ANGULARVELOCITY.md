# Velocidad Angular (Omega) desde preCICE

## Descripción

El módulo FSI ahora soporta recibir la **velocidad angular (Omega)** desde preCICE y almacenarla en un campo escalar `UniformDimensionedScalarField` de OpenFOAM. Este campo puede ser utilizado por otros modelos, condiciones de borde, y funciones que necesiten la velocidad de rotación.

## Motivación

En OpenFOAM v2406+, se introdujo la capacidad de usar campos `UniformDimensionedScalarField` con la función `lookup` en condiciones de borde y modelos de movimiento. Esto permite que objetos rotativos utilicen una velocidad angular calculada dinámicamente, por ejemplo desde un solver estructural acoplado.

## ¿Cómo Funciona?

### Flujo de Datos

```
preCICE (Solver Estructural/Controlador)
    ↓
AngularVelocity [scalar] → Omega en rad/s
    ↓
Interface.readCouplingData()
    ↓
AngularVelocity.read()
    ↓
uniformDimensionedScalarField "omega"
    ↓
Disponible para:
  - solidBodyMotion (rotatingMotion)
  - Boundary conditions
  - Function objects
  - Custom models
```

### Campo Creado

El adaptador crea/actualiza un `uniformDimensionedScalarField` con:
- **Nombre**: Configurable (por defecto `"omega"`)
- **Dimensiones**: `[0 0 -1 0 0 0 0]` (rad/s)
- **Ubicación**: Registrado en el object registry del mesh
- **Accesibilidad**: Global, disponible mediante `lookup`

## Configuración

### 1. En preciceDict

```cpp
FSI
{
    solverType incompressible;
    
    // Nombre del campo de velocidad angular (opcional, por defecto "omega")
    nameOmegaField omega;
    
    // Eje de rotación (para referencia, si se necesita vector omega)
    rotationAxis (0 0 1);  // Eje Z
}
```

### 2. En precice-config.xml

```xml
<participant name="Fluid">
    <use-mesh name="Fluid-Mesh" provide="yes"/>
    
    <!-- Escribir fuerzas/torque -->
    <write-data name="Force" mesh="Fluid-Mesh"/>
    <write-data name="Torque" mesh="Fluid-Mesh"/>
    
    <!-- Leer velocidad angular -->
    <read-data name="AngularVelocity" mesh="Fluid-Mesh"/>
</participant>

<participant name="Structure">
    <use-mesh name="Fluid-Mesh" from="Fluid"/>
    
    <!-- Recibir fuerzas/torque -->
    <read-data name="Force" mesh="Fluid-Mesh"/>
    <read-data name="Torque" mesh="Fluid-Mesh"/>
    
    <!-- Enviar velocidad angular calculada -->
    <write-data name="AngularVelocity" mesh="Fluid-Mesh"/>
</participant>
```

### 3. En system/preciceDict

```cpp
preciceConfig "precice-config.xml";
participant Fluid;
modules (FSI);

interfaces
{
    RotorInterface
    {
        mesh Fluid-Mesh;
        patches (rotor);
        locations faceCenters;
        
        writeData
        (
            Force
            Torque
        );
        
        readData
        (
            AngularVelocity
        );
    };
};

FSI
{
    solverType incompressible;
    
    // Campo para almacenar omega
    nameOmegaField omega;
    
    // Eje de rotación
    rotationAxis (0 0 1);
    
    // Densidad y viscosidad
    rho rho [1 -3 0 0 0 0 0] 1000;
    nu nu [0 2 -1 0 0 0 0] 1e-6;
}
```

## Uso del Campo Omega en OpenFOAM

### Ejemplo 1: solidBodyMotion con rotatingMotion

En `constant/dynamicMeshDict`:

```cpp
dynamicFvMesh   solidBodyMotionFvMesh;

motionSolverLibs ("libfvMotionSolvers.so");

solidBodyMotionFunction rotatingMotion;

rotatingMotionCoeffs
{
    origin      (0 0 0);
    axis        (0 0 1);
    
    // Usar omega desde preCICE
    omega       lookup;
    name        omega;  // Nombre del campo (debe coincidir con nameOmegaField)
}
```

### Ejemplo 2: Condición de Borde con Velocidad Variable

En `0/U`:

```cpp
rotatingWall
{
    type            rotatingWallVelocity;
    origin          (0 0 0);
    axis            (0 0 1);
    
    // Usar omega desde preCICE
    omega           lookup;
    name            omega;
}
```

### Ejemplo 3: Function Object que Usa Omega

En `system/controlDict`:

```cpp
functions
{
    tipSpeed
    {
        type            coded;
        name            calculateTipSpeed;
        
        codeExecute
        #{
            // Leer omega desde el field registry
            const auto& omega = mesh().lookupObject<uniformDimensionedScalarField>("omega");
            
            scalar radius = 0.5;  // Radio del rotor
            scalar tipSpeed = omega.value() * radius;
            
            Info << "Tip speed: " << tipSpeed << " m/s" << endl;
            Info << "Omega: " << omega.value() << " rad/s ("
                 << omega.value() * 30.0 / Foam::constant::mathematical::pi << " RPM)" << endl;
        #};
    }
}
```

## Formato de Datos preCICE

### AngularVelocity

- **Tipo**: Escalar
- **Unidad**: rad/s
- **Valores**: Una velocidad angular por cada punto de la interfaz (típicamente el mismo valor)

Ejemplo desde el participante estructural (Python):

```python
import numpy as np

# Calcular velocidad angular del sistema estructural
# Ejemplo: de un controlador PID basado en torque
torque_from_fluid = precice.read_data(...)
omega_rad_s = calculate_angular_velocity(torque_from_fluid)

# Convertir a RPM si es útil para debugging
omega_rpm = omega_rad_s * 30.0 / np.pi
print(f"Omega: {omega_rad_s:.4f} rad/s ({omega_rpm:.2f} RPM)")

# Crear array con el mismo valor para todos los puntos
num_vertices = precice.get_mesh_vertex_size(mesh_id)
omega_data = np.full(num_vertices, omega_rad_s)

# Escribir a preCICE
precice.write_data(mesh_id, "AngularVelocity", vertex_ids, omega_data)
```

## Casos de Uso

### 1. Turbina con Control de Velocidad

**Escenario**: Una turbina eólica donde el controlador ajusta la velocidad de rotación basándose en las condiciones del viento.

```
Fluido → Calcula fuerzas aerodinámicas
    ↓
Estructura → Calcula torque neto
    ↓
Controlador → Ajusta omega basándose en torque
    ↓
Fluido → Recibe nueva omega para rotación
```

**Configuración**:
```cpp
readData (AngularVelocity);
writeData (Force Torque);

// En dynamicMeshDict
omega lookup;
name omega;
```

### 2. Motor con Acoplamiento Termomecánico

**Escenario**: Un motor donde la velocidad depende de la temperatura y carga.

```
Fluido → Calcula flujo térmico y fuerzas
    ↓
Estructura → Calcula expansión térmica y respuesta mecánica
    ↓
Modelo de motor → Calcula omega(temperatura, carga)
    ↓
Fluido → Ajusta rotación según omega
```

### 3. Bomba con Control PID

**Escenario**: Bomba centrífuga con controlador que mantiene presión objetivo.

```
Fluido → Mide presión de salida
    ↓
Controlador PID → Ajusta omega para mantener presión
    ↓
Fluido → Aplica nueva velocidad de rotación
```

## Ventajas

✅ **Compatible**: Usa el sistema nativo de OpenFOAM v2406+ (`lookup`)  
✅ **Flexible**: El campo está disponible globalmente  
✅ **Eficiente**: Campo escalar uniforme (memoria mínima)  
✅ **Reutilizable**: Múltiples modelos pueden usar el mismo campo  
✅ **Debugging fácil**: El campo es visible y puede ser monitoreado  

## Diferencias con RotationAngle

| Aspecto | AngularVelocity | RotationAngle |
|---------|-----------------|---------------|
| **Tipo** | Velocidad (rad/s) | Posición angular (rad) |
| **Uso** | Para modelos que necesitan ω | Para aplicar rotación directa |
| **Almacenamiento** | UniformDimensionedScalarField | No almacena (aplica directamente) |
| **Aplicación** | Usado por OpenFOAM internamente | Transforma mesh_.points() |
| **Compatibilidad** | solidBodyMotion, BC, etc. | Rotación directa de malla |

## Combinación con Otros Datos

Puedes combinar `AngularVelocity` con otros datos FSI:

```cpp
readData
(
    Displacement       // Deformación flexible
    RotationAngle      // Posición angular absoluta
    AngularVelocity    // Velocidad de rotación
);
```

**Uso típico combinado:**
1. `RotationAngle` → Rota toda la malla a la posición actual
2. `Displacement` → Agrega deformación elástica
3. `AngularVelocity` → OpenFOAM usa para cálculos internos (BC, fuerzas centrífugas, etc.)

## Ejemplo Completo: Turbina con Control

### Configuración preCICE

```cpp
FSI
{
    solverType incompressible;
    
    // Velocidad angular
    nameOmegaField omega;
    rotationAxis (0 0 1);
    
    // Rotación
    rotationCenter (0 0 0);
    
    // Torque
    torqueReferencePoint (0 0 0);
}

interfaces
{
    TurbineRotor
    {
        mesh Rotor-Mesh;
        patches (blades);
        
        writeData
        (
            Force
            Torque
        );
        
        readData
        (
            RotationAngle     // Posición angular
            AngularVelocity   // Velocidad de rotación
        );
    };
};
```

### dynamicMeshDict

```cpp
dynamicFvMesh solidBodyMotionFvMesh;
motionSolverLibs ("libfvMotionSolvers.so");

solidBodyMotionFunction rotatingMotion;

rotatingMotionCoeffs
{
    origin (0 0 0);
    axis (0 0 1);
    omega lookup;
    name omega;  // Campo desde preCICE
}
```

### Resultado

- Las fuerzas aerodinámicas se calculan correctamente
- El torque se envía al controlador
- El controlador ajusta omega
- OpenFOAM usa la nueva omega para:
  - Rotación de la malla
  - Cálculo de fuerzas centrífugas
  - Condiciones de borde rotatorias

## Debugging

### Verificar el Campo en Runtime

```cpp
// En cualquier código (function object, BC, etc.)
const auto& omega = mesh.lookupObject<uniformDimensionedScalarField>("omega");
Info << "Current omega: " << omega.value() << " rad/s" << endl;
```

### En el Log del Adaptador

Con `DEBUG` activo, verás:

```
AngularVelocity coupling initialized:
  Field name: omega
  Rotation axis: (0, 0, 1)
Received angular velocity: 10.472 rad/s (100 RPM)
Updated field 'omega' = 10.472 rad/s
```

### Verificar Disponibilidad

```cpp
if (mesh.foundObject<uniformDimensionedScalarField>("omega"))
{
    Info << "Omega field is available" << endl;
}
else
{
    Info << "Omega field NOT found" << endl;
}
```

## Notas Importantes

⚠️ **Unidades**: Siempre en rad/s, no en RPM

⚠️ **OpenFOAM versión**: La función `lookup` requiere OpenFOAM v2406 o posterior

⚠️ **Nombre del campo**: Debe coincidir entre `nameOmegaField` y el uso en OpenFOAM

⚠️ **Registro**: El campo se registra automáticamente en el mesh object registry

⚠️ **Persistencia**: El campo existe mientras el caso está corriendo, no se guarda en disco por defecto

## Conversión de Unidades

```cpp
// rad/s a RPM
RPM = omega * 30.0 / π

// RPM a rad/s
omega = RPM * π / 30.0

// Ejemplo: 100 RPM
omega = 100 * 3.14159 / 30.0 = 10.472 rad/s
```
