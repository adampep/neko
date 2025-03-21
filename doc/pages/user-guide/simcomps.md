# Simulation components {#simcomps}

\tableofcontents

## What are simulation components?
Simulation components, or simcomps fo short,  incapsulate additional
functionality that may be useful for certain cases but not necessary to run the
solver.
This can include computation and output of additional fields, in-situ
post-processing operations, data sampling, etc.

By design, simulation components can tap into every aspect of the simulation,
so they can be quite powerful.
As the code grows, we expect to add more and more simcomps to the code.

## Adding simulation components to the case
Each simcomp is defined as a single JSON object at are added to an array
of objects called `simulation_components`, which resides directly under the
`case` object.

## List of simulation components

The following is a list of simulation components that are currently available
in Neko. The list will be updated as new simcomps are added.

- Computation of vorticity field \ref simcomp_vorticity
- Computation of \f$ \lambda_2 \f$ \ref simcomp_lambda2
- Probing of fields at selected points \ref simcomp_probes
- Output of registered fields to an `.fld` file \ref simcomp_field_writer
- Computation of the derivative of a field \ref simcomp_derivative
- Computation of forces and torque on a surface \ref simcomp_force_torque
- Computation of the weak gradient of a field \ref simcomp_weak_grad
- Computation of subgrid-scale (SGS) eddy viscosity via a SGS model \ref simcomp_les_model
- User defined components \ref user-file_simcomps
- Fluid statistics simcomp, "fluid_stats", for more details see the [statistics guide](@ref statistics-guide)
- Computation of the spectral error indicator \ref simcomp_speri

## Controling execution and file output
Each simulation component is, by default, executed once per time step to perform
associated computations and output. However, this can be modified by using the
`compute_control` and `compute_value` parameters for the computation and the
`output_control` and `output_value` for the output to disk. The parameters
for the `_control` values are the same as for the fluid and checkpointing.
Additionally, one can set `output_control` to `global` and `never`. The former
will sync the `output_` parameter to that of the fluid. Choosing `never` will
suppress output all together. If no parameters for the `output_` parameters are
provided, they are set to be the same as for `compute_`. In order to simplify
the configuration, the `compute_control` can be set to `fluid_output` to sync
the computation to the fluid output.

For simcomps that compute 3D fields, the output can be either added to the main
`.fld` file, containing velocity and pressure, or saved to a separate file. For
the latter, the `output_filename` keyword should be provided. One can
additionally provide the `precision` keyword, which can be set to either
`single` or `double` to control the precision of the written data.

For example, in the `tgv` example case the `vorticity` component is executed
once per 50 time steps. The `output_` parameters are synced to that, and the
vorticity fields will be added to the main `.fld` file.
~~~~~~~~~~~~~~~{.json}
{
    "type": "vorticity",
    "compute_control": "tsteps",
    "compute_value": 50
}
~~~~~~~~~~~~~~~

### vorticity {#simcomp_vorticity}
Computes the vorticity field an stores in the field registry as `omega_x`,
`omega_y` and `omega_z`. By default, appends the 3 vorticity fields to the field files as
scalars. To output in a different `fld` series, use the `"output_filename"` parameter.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "vorticity"
 }
 ~~~~~~~~~~~~~~~

### lambda2 {#simcomp_lambda2}
Computes \f$ \lambda_2 \f$ for the velocity field and stores it in the normal output files as the first unused field.
This means that \f$ \lambda_2 \f$ can be found in the temperature field in then fld files if running without a scalar
and s1 if neko is run with one scalar. To output in a different `fld` series, use the `"output_filename"` parameter.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "lambda2"
 }
 ~~~~~~~~~~~~~~~

### probes {#simcomp_probes}
Probes selected solution fields at a list of points. This list of points can be
generated in a variety of ways, but the most common is to use the `csv` type.

#### Supported types:

 - `file`: Reads a list of points from a CSV file. The name of the file is
   provided with the `file_name` keyword. The CSV file should have the
   following format:
   ~~~~~~~~~~~~~~~{.csv}
   x_0, y_0, z_0
   x_1, y_1, z_1
   ...
   x_N, y_N, z_N
   ~~~~~~~~~~~~~~~
   The points are assumed to be in the same units as the simulation.
- `points`: Reads a list of points from a JSON file. The points are specified
  based in the `coordinates` keyword and should be a list of x,y,z values.
  The file should have the following format:
  ~~~~~~~~~~~~~~~{.json}
  {
    "type": "points",
    "coordinates": [
      0.0, 0.0, 0.0,
      1.0, 0.0, 0.0,
      0.0, 1.0, 0.0,
      0.0, 0.0, 1.0,
      ...
    ]
  }
  ~~~~~~~~~~~~~~~
  The points are assumed to be in the same units as the simulation.
- `line`: Generates a list of points along a line. The line is defined by two
  points, `start` and `end`, and the number of points to generate, `amount`.
  The points are generated by linearly interpolating between `start` and `end`.
  The line is defined as:
  ~~~~~~~~~~~~~~~{.json}
  {
    "type": "line",
    "start": [0.0, 0.0, 0.0],
    "end": [1.0, 1.0, 1.0],
    "amount": 10
  }
  ~~~~~~~~~~~~~~~
- `circle`: Generates a list of points along a circle. The circle is defined by
  a center, `center`, a radius, `radius` and the normal, `normal`, the number of
  points to generate is controlled by `amount`. The points are generated by
  rotating a point around the center starting from the specified axis projected
  onto the circle. The circle is defined:
  ~~~~~~~~~~~~~~~{.json}
  {
    "type": "circle",
    "center": [0.0, 0.0, 0.0],
    "radius": 1.0,
    "normal": [0.0, 0.0, 1.0],
    "axis": "x",
    "amount": 4
  }
  ~~~~~~~~~~~~~~~
  Leads to the following points:
  ~~~~~~~~~~~~~~~{.csv}
  1.0, 0.0, 0.0
  0.0, 1.0, 0.0
  -1.0, 0.0, 0.0
  0.0, -1.0, 0.0
  ~~~~~~~~~~~~~~~

 #### Example usage:
 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "probes",
   "compute_control": "simulationtime",
   "compute_value"    : 1,
   "fields": ["w","s"],
   "output_file":  "output.csv",
   "points": [
      {
        "type": "file",
        "file_name": "points.csv"
      }
    ],
 }
 ~~~~~~~~~~~~~~~
This probes the fields 'w', and 's' in the points described by points.csv and
outputs into output.csv every 1 time units.

The probed information will be saved in the output file in the following format:

~~~~~~~~~~~~~~~{.csv}
N_p, N_f, fields[0], fields[1], ..., fields[N_f-1]
p_0_x, p_0_y, p_0_z
p_1_x, p_1_y, p_1_z
...
p_N_p_x, p_N_p_y, p_N_p_z
time_0, p_0_field_0, p_0_field_1, ..., p_0_field_N_f-1
time_0, p_1_field_0, p_1_field_1, ..., p_1_field_N_f-1
...
time_0, p_N_p_field_0, p_N_p_field_1, ..., p_N_p_field_N_f-1
time_1, p_0_field_0, p_0_field_1, ..., p_0_field_N_f-1
time_1, p_1_field_0, p_1_field_1, ..., p_1_field_N_f-1
...
time_N_p, p_N_p_field_0, p_N_p_field_1, ..., p_N_p_field_N_f-1
~~~~~~~~~~~~~~~

### field_writer {#simcomp_field_writer}
Outputs registered 3D fields to an `.fld` file. Requires a list of field names
in the `fields` keyword. Primarily to be used for outputting new fields defined
in the user file. The fields are added to then `neko_field_registry` object and
are expected to be updated in the user file, or, perhaps, by other simcomps.
Since this simcomp does not compute anything `compute_` configuration is
irrelevant.
 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "field_writer",
   "fields": ["my_field1", "my_field2"],
   "output_filename": "myfields",
   "precision": "double",
   "output_control" : "simulation_time",
   "output_value" : 1.0
 }
 ~~~~~~~~~~~~~~~

### derivative {#simcomp_derivative}
Computes the derivative of field along a chosen direction (x, y, or z). The
field to derivate is controlled by the `field` keyword and the direction by the
`direction` keyword. The simcomp will register the computed derivatives in the
registry as `d[field]_d[direction]`, where the values in the brackets
correspond to the choice of the user keywords. Supports writing the computed
fields to disk via the usual common keywords. The resulting field will be
appended as a scalar to the field files. To output in a different `fld` series,
use the `"output_filename"` parameter.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "derivative",
   "field": "u",
   "direction": "y"
 }
 ~~~~~~~~~~~~~~~

### force_torque {#simcomp_force_torque}
Computes the force on a specified zone and the corresponding torque
around a center point. The compute control specifies how often they are
computed and printed into the log. Scale specifies a scale for the computed
force/torque. Conventient if one wants to scale with the area or similar. long_print is default false and can be set to true to print all digits in the calculation.
Subroutines used in the simcomp can be found in src/qoi/drag_torque.f90

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "force_torque",
   "zone_id": 1,
   "center": [0.0, 0.0, 0.0],
   "zone_name": "some chosen name, optional",
   "scale": 1.0
   "long_print" : false
   "compute_control" : "tsteps",
   "compute_value" : 10
 }
 ~~~~~~~~~~~~~~~

### weak_grad {#simcomp_weak_grad}
Computes the weak gradient of a field. The weak gradient is value of the
gradient multiplied by the local value of the mass matrix. This is how a
gradient term appears in the weak formulation of the governing equations. The
field to derivate is controlled by the `field` keyword. The simcomp will
register the computed components of the gradients in the registry as
`weak_grad_[field]_x`, `weak_grad_[field]_y`, `weak_grad_[field]_z` where the
value in the brackets corresponds to the choice of the user keyword. Supports
writing the computed fields to disk via the usual common keywords.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "weak_grad"
   "field": "u",
   "output_control" : "never"
 }
 ~~~~~~~~~~~~~~~

### les_model {#simcomp_les_model}
Computes a subgrid eddy viscosity field using an SGS model. **Note*:* The simcomp
*only* computes the eddy viscosity field. You have to select the corresponding
`nut_field` in the fluid and/or scalar JSON object to actually enable LES, see
corresponding documentation. The simcomp is controlled by the following
keywords:

- `model`: Selects the SGS model. Currently available models are:
  - `smagorinsky`: The standard Smagorinsky model. Configured by the
    following additional keyword:
    - `c_s`: The Smagorinsky constant, defaults to 0.17.
  - `dynamic_smagorinsky`: The dynamic Smagorinsky model.
  - `vreman`: The Vreman model. Configured by the following additional keyword:
    - `c`: The model constant, defaults to 0.07.
  - `sigma`: The Sigma model. Configured by the following additional keyword:
    - `c`: The model constant, defaults to 1.35.
  - `wale`: The WALE model. Configured by the following additional keyword:
    - `c_w`: The WALE constant, defaults to 0.55.
- `les_delta`: Selects the way to compute the LES filter length scale. Currently three
  alternatives are provided and the default one is `pointwise` if
  nothing is specified:
  - `pointwise`: Computes a local value based on the spacing of the GLL nodes.
  - `elementwise_average`: Computes a single value for the whole element based on the
    average spacing of the GLL nodes within the element.
  - `elementwise_max`: Computes a single value for the whole element based on the
    maximum spacing of the GLL nodes within the element.
  The `les_delta` field is added to the registry and written to the .fld files.
- `nut_field`: The name of the SGS eddy viscosity field added to the registry.
  Defaults to `nut`. This allows to have two different SGS models active, saved
  to different fields. For example, one for the scalar and one to the fluid.
- `extrapolation`: Whether or not extrapolate the velocity to
  compute the eddy viscosity.
  - `true`: extrapolate the velocity as the same order as
  the time scheme.
  - `false`: the default option, disable the extrapolation. 
  In this case, the estimation of the eddy viscosity is of first order, while 
  circumvent the risk of unstable extrapolation.

 ~~~~~~~~~~~~~~~{.json}
 {
   "type": "les_model"
   "model": "smagorinsky",
   "delta_type": "pointwise",
   "output_control" : "never"
 }
 ~~~~~~~~~~~~~~~

### Spectral error indicator {#simcomp_speri}

Computes the spectral error indicator as developed by Mavriplis (1989) (https://doi.org/10.1007/978-3-663-13975-1_34).
This is an a posteriori error measure, based on the local properties of
the spectral solution. This method formally only gives an indication of the error.

The spectral error indicator is computed for the 3 velocity fields, resulting
in 3 additional fields appended to the field files.

~~~~~~~~~~~~~~~{.json}
 {
   "type": "spectral_error"
 }
 ~~~~~~~~~~~~~~~
