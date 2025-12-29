set DEBRIS; # Set of debris objects

param N := card(DEBRIS);  # Number of debris objects

param orbit {DEBRIS} integer >= 0;  # Orbit altitude (km)

param m0 integer >= 600;      # Initial mothership mass (kg) (600 kg is the minimum mass for the mothership)
param m_n integer >= 0;       # Count of de-orbiting spacecraft
param m_sc integer >= 0;      # Mass of de-orbiting spacecraft (kg)
param Isp integer >= 0;       # Specific impulse (s)
param fuel0 integer >= 0;     # Initial available fuel (kg)
param mu integer := 3.986e14; # Gravitational constant (m^3/s^2)
param g0 := 9.81;             # Standard gravity (m/s^2)
param o0 integer >= 0;        # Initial mothership orbit (km)

param value {DEBRIS} integer >= 1, <= 10; # Value of each debris object (arbitrary units)

param transfer_DeltaV {i in DEBRIS, j in DEBRIS: i != j} := # Delta-V required to transfer from debris i orbit to debris j orbit (m/s) 
    if i = j then 0
    else  
        sqrt(mu / orbit[i] * 1.0e3) * (sqrt(2 * orbit[j] * 1.0e3 / (orbit[i] * 1.0e3 + orbit[j] * 1.0e3)) - 1) + 
        sqrt(mu / orbit[j] * 1.0e3) * (1 - sqrt(2 * orbit[i] * 1.0e3 / (orbit[i] * 1.0e3 + orbit[j] * 1.0e3))); 

var x {i in DEBRIS, j in DEBRIS} binary; # 1 if debris i is transferred to orbit of debris j, 0 otherwise
var y {i in DEBRIS} binary; # 1 if debris i is de-orbited, 0 otherwise
var u {i in DEBRIS} integer >= 1, <= N; # Auxiliary variable to ensure each debris is visited only once

var mass {k in 0..N};      # Mass of mothership after k de-orbiting spacecraft removal (kg)
var fuel_used {k in 0..N}; # Total fuel used (kg)

maximize total_value: 
    sum {i in DEBRIS} value[i] * y[i]; # Maximize total value of de-orbited debris

subject to Initial_mass:
    mass[0] = m0; # Initial mass of mothership

subject to Fuel_calculation{k in 1..N}:
    fuel_used[k] = sum {i in DEBRIS, j in DEBRIS: i != j} (0.5 * mass[k - 1] * exp(-(transfer_DeltaV[i,j]/(Isp * g0)))); # Total fuel used for transfers 

subject to Mass_Update {k in 1..N}:
    mass[k] = mass[k-1] - m_sc - fuel_used[k]; # Update mass after each de-orbiting spacecraft removal and transfer

var is_start {i in DEBRIS} binary;
var is_end {i in DEBRIS} binary;

subject to One_Start:
    sum {i in DEBRIS} is_start[i] = 1;

subject to One_End:
    sum {i in DEBRIS} is_end[i] = 1;

subject to Incoming_Arc {j in DEBRIS}:
    sum {i in DEBRIS: i != j} x[i,j] = y[j] - is_start[j];

subject to Outgoing_Arc {i in DEBRIS}:
    sum {j in DEBRIS: j != i} x[i,j] = y[i] - is_end[i];

subject to Subtour_Elimination {i in DEBRIS, j in DEBRIS: i != j}:
    u[i] + 1 <= u[j] + N * (1 - x[i,j]);

subject to Remove_count:
    sum {i in DEBRIS} y[i] <= m_n; # Total number of de-orbiting spacecraft must be less or equal to m_n

subject to Fuel_Constraint:
    sum {i in 1..N} fuel_used[i] <= fuel0; # Total fuel used must not exceed initial available fuel


