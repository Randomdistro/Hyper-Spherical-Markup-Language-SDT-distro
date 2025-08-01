# HSML Language Specification v21.0

## Overview

HSML (Hyper-Spherical Markup Language) is a markup language designed for describing 3D scenes using pure spherical coordinates and Spatial Displacement Theory (SDT) physics. It eliminates all Cartesian contamination and operates exclusively in spherical truth.

## Core Philosophy

1. **Pure Spherical Coordinates**: All positioning uses (r, θ, φ) with zero-division safe 1-361° range
2. **21-Dimensional State**: Every element maintains a 21D state vector
3. **No Forces**: Only spatial displacement according to SDT principles
4. **Resonance-Based**: All interactions based on harmonic resonance patterns

## Document Structure

```hsml
<?hsml version="21.0" encoding="spherical-utf8"?>
<hsml:sphere xmlns:hsml="https://hsml.sdt/2024" sdt="authentic">
  <!-- Document content -->
</hsml:sphere>
```

## Root Elements

### `<hsml:sphere>`
The root container representing the spherical universe.

**Attributes:**
- `sdt` (required): SDT implementation level (`"authentic"`, `"compatible"`, `"emulated"`)
- `resonance`: Base resonance frequency in Hz (default: `432`)
- `time-scale`: Time dilation factor (default: `1.0`)

```hsml
<hsml:sphere sdt="authentic" resonance="528" time-scale="1.2">
  <!-- Content -->
</hsml:sphere>
```

### `<hsml:viewport>`
Defines the steradian viewport with four-corner interpolation.

**Attributes:**
- `steradian-corners`: Number of interpolation corners (must be `4`)
- `interpolation`: Interpolation method (`"spherical-native"`, `"linear"`)
- `bubble-radius`: Viewport bubble radius (default: `∞`)
- `follow-user`: Whether viewport follows user position (`true`/`false`)

```hsml
<hsml:viewport 
    steradian-corners="4"
    interpolation="spherical-native"
    bubble-radius="100.0"
    follow-user="true">
  <!-- Scene content -->
</hsml:viewport>
```

## Positioning Elements

### `<hsml:element>`
Basic positioned element in spherical space.

**Position Attributes:**
- `r`: Radius from origin (positive real number)
- `θ` (theta): Polar angle in degrees (1-361, safe range)
- `φ` (phi): Azimuthal angle in degrees (1-361, safe range)

**Physics Attributes:**
- `state`: Matter state (`"solid"`, `"liquid"`, `"gas"`, `"plasma"`, `"quantum"`, `"void"`)
- `resonance`: Element resonance frequency in Hz
- `mass-equiv`: Mass equivalent for SDT calculations
- `collision-radius`: Spherical collision boundary

**21D State Attributes:**
- `zero-point`: Level 0 dimension (1D)
- `line-x`, `line-y`: Level 1 dimensions (2D)
- `plane-x`, `plane-y`, `plane-z`: Level 2 dimensions (3D)
- `sphere-r`, `sphere-θ`, `sphere-φ`, `sphere-t`: Level 3 dimensions (4D)
- `flux-*`: Level 4 dimensions (5D)
- `energy-*`: Level 5 dimensions (6D)

```hsml
<hsml:element 
    r="10.5" θ="180" φ="90"
    state="plasma"
    resonance="432"
    mass-equiv="1.0"
    collision-radius="1.0"
    zero-point="1.0"
    sphere-r="10.5">
  <!-- Element content -->
</hsml:element>
```

## Specialized Elements

### `<hsml:field>`
Defines SDT fields that influence spatial displacement.

**Attributes:**
- `type`: Field type (`"displacement"`, `"resonance"`, `"flux"`, `"temporal"`)
- `strength`: Field strength multiplier
- `radius`: Field influence radius
- `decay`: Field decay function (`"linear"`, `"quadratic"`, `"exponential"`)

```hsml
<hsml:field 
    type="displacement"
    r="0" θ="180" φ="180"
    strength="5.0"
    radius="50.0"
    decay="quadratic">
</hsml:field>
```

### `<hsml:user>`
Represents the user as a physics object for collision detection.

**Attributes:**
- `collision`: Enable user collision (`true`/`false`)
- `matter-state`: User's matter state (affects interactions)
- `resonance`: User's resonance frequency

```hsml
<hsml:user 
    r="5.0" θ="180" φ="90"
    collision="true"
    matter-state="quantum"
    resonance="528">
</hsml:user>
```

### `<hsml:resonator>`
Creates resonance patterns affecting nearby elements.

**Attributes:**
- `frequency`: Primary resonance frequency
- `harmonics`: Comma-separated harmonic frequencies
- `pattern`: Resonance pattern (`"sine"`, `"square"`, `"sawtooth"`, `"rrpt"`)
- `amplitude`: Resonance amplitude

```hsml
<hsml:resonator
    r="15.0" θ="225" φ="135"
    frequency="432"
    harmonics="864,1296"
    pattern="rrpt"
    amplitude="2.0">
</hsml:resonator>
```

## Content Elements

### Text Content
Text is rendered in 3D space using spherical typography.

```hsml
<hsml:text r="8.0" θ="180" φ="90" size="2.0" color="#FF6B35">
  Spherical Truth Revealed
</hsml:text>
```

### Media Elements

```hsml
<hsml:model 
    r="12.0" θ="200" φ="75"
    src="assets/sphere.obj"
    format="spherical-obj">
</hsml:model>

<hsml:texture
    r="5.0" θ="160" φ="120"
    src="assets/pattern.png"
    projection="spherical">
</hsml:texture>
```

## Interaction Elements

### `<hsml:trigger>`
Defines spatial triggers that respond to user presence.

**Attributes:**
- `shape`: Trigger shape (`"sphere"`, `"steradian-cone"`)
- `size`: Trigger size
- `event`: Event type (`"enter"`, `"exit"`, `"dwell"`)

```hsml
<hsml:trigger
    r="20.0" θ="180" φ="90"
    shape="sphere"
    size="3.0"
    event="enter"
    action="displace-to(r=25, θ=180, φ=90)">
</hsml:trigger>
```

### `<hsml:animation>`
Defines spatial displacement animations.

```hsml
<hsml:animation
    target="#element1"
    type="displacement"
    duration="5.0"
    easing="spherical-smooth">
  
  <hsml:keyframe time="0.0" r="10" θ="180" φ="90"/>
  <hsml:keyframe time="1.0" r="15" θ="225" φ="135"/>
  
</hsml:animation>
```

## Styling

### HSML Stylesheets (.hss files)

```hss
/* Spherical Stylesheet */
@sphere-rule {
    element[state="plasma"] {
        resonance: 528Hz;
        flux-density: 0.8;
        matter-coherence: quantum;
        glow-intensity: 2.0;
    }
    
    field[type="displacement"] {
        visualization: steradian-field-lines;
        color: rgba(255, 107, 53, 0.3);
    }
    
    viewport {
        steradian-coverage: 4π;
        interpolation-quality: authentic;
        background: void-black;
    }
}
```

## Scripting

### Spherical Script (.ss files)

```ss
// Spherical Script
namespace MyScene {
    sphere main() {
        let origin = SphericalCoord(r: 1.0, θ: 180, φ: 90);
        let user_pos = get_user_position();
        
        on_user_move(position: SphericalCoord) {
            if (distance(position, origin) < 5.0) {
                emit_resonance(frequency: 432Hz);
            }
        }
        
        function calculate_displacement(from: SphericalCoord, to: SphericalCoord) -> SphericalCoord {
            return SphericalCoord(
                r: (to.r - from.r) * 0.1,
                θ: safe_angle(to.θ - from.θ),
                φ: safe_angle(to.φ - from.φ)
            );
        }
    }
}
```

## Comments and Documentation

```hsml
<!-- Standard HSML comment -->
<hsml:comment type="physics">
    This element exists in pure spherical space
    at radius 10.5, polar angle 180°, azimuthal angle 90°
    with authentic SDT physics enabled
</hsml:comment>
```

## Validation and Safety

### Zero-Division Safety
All angles are automatically clamped to the safe 1-361° range:

```hsml
<!-- These are automatically made safe -->
<hsml:element r="5.0" θ="0" φ="365"/>
<!-- Becomes: θ="1" φ="361" -->
```

### Physics Validation
Elements with invalid 21D states are rejected:

```hsml
<!-- Invalid: infinite energy -->
<hsml:element r="5.0" energy-quantum="∞"/>
<!-- Error: Invalid 21D state -->
```

### SDT Compliance Levels

1. **`sdt="authentic"`**: Full SDT physics, no Cartesian conversion
2. **`sdt="compatible"`**: SDT physics with Cartesian fallback for legacy
3. **`sdt="emulated"`**: Cartesian simulation of SDT behavior

## File Extensions

- `.hsml`: HSML markup files
- `.hss`: HSML stylesheet files
- `.ss`: Spherical Script files
- `.hsml-project`: Project configuration files

## MIME Types

- `application/hsml+xml`: HSML markup
- `text/hsml-stylesheet`: HSML stylesheets
- `application/spherical-script`: Spherical Script

## Example Complete Document

```hsml
<?hsml version="21.0" encoding="spherical-utf8"?>
<hsml:sphere xmlns:hsml="https://hsml.sdt/2024" sdt="authentic" resonance="432">
  
  <hsml:viewport 
      steradian-corners="4"
      interpolation="spherical-native"
      bubble-radius="∞"
      follow-user="true">
    
    <!-- Central resonator -->
    <hsml:resonator
        r="0" θ="180" φ="180"
        frequency="432"
        pattern="rrpt"
        amplitude="1.0">
    </hsml:resonator>
    
    <!-- Orbiting plasma element -->
    <hsml:element 
        id="plasma-orb"
        r="10.0" θ="180" φ="90"
        state="plasma"
        resonance="864"
        collision-radius="1.5">
      
      <hsml:text size="1.0" color="#FF6B35">
        Spherical Truth
      </hsml:text>
      
    </hsml:element>
    
    <!-- Displacement field -->
    <hsml:field
        type="displacement"
        r="0" θ="180" φ="180"
        strength="2.0"
        radius="20.0"
        decay="quadratic">
    </hsml:field>
    
    <!-- User representation -->
    <hsml:user
        collision="true"
        matter-state="quantum"
        resonance="528">
    </hsml:user>
    
  </hsml:viewport>
  
</hsml:sphere>
```

This specification defines a complete markup language for spherical 3D content with authentic SDT physics, eliminating all Cartesian contamination and providing a pure spherical development environment.