# ShW-MAGICA
ShW-MAGICA is a Python implementation of the Lehtinen–Pirjola method for computing Geomagnetically Induced Currents (GICs) in transmission networks, including the explicit modelling of shield wires (ShW).

The code is based on the original GEOMAGICA framework and extends it by incorporating the LPm method, explicit shield wire equivalent circuits, and adaptations for the Portuguese transmission network.

# Input Files

GRID.txt - Defines the network nodes. Each row corresponds to one node.
Column	Description: Node number, Node name, Node code,	Country,	Latitude (°),	Longitude (°), Grounding resistance (Ω), Transformer resistance (Ω)

CONNECTIONS.txt - Defines the network connections. Each row corresponds to one line or transformer.
Column	Description: substation name, transformer type, winding, from node, to node, line length, shield wire type, phase exception, line resistance, nominal voltage

GRID_ShW.txt - Defines the network nodes regarding ShW. Each row corresponds to one new node.
Column	Description: node identifier, latitude, longitude, grounding resistance, transformer resistance

CONNECTIONS_ShW.txt - Defines the network connections regarding ShW. Each row corresponds to one new line of ShW.
Column	Description: shield wire node, connected power system node, corresponding transmission line, equivalent resistance, grounding parameters

Point_Coord.pkl - Contains the coordinates of the electric field interpolation grid
The file stores two arrays: latitude, longitude

E_field_Coord.pkl - The electric field file contains the geoelectric field evaluated at the interpolation points defined in Point_Coord.pkl.
Each file stores both field components: northward component (E_n), eastward component (E_e)

# Credits

This work is based on the original GEOMAGICA framework developed from

P. Weidelt, C. Beggan, K. Turnbull, A. McKay, R. Bailey
https://github.com/geomagpy/GEOMAGICA

The present implementation (ShW-MAGICA) was developed by

Rute Rodrigues dos Santos, Alexandra Pais, Joana Alves Ribeiro, Fernando Pinheiro

Department of Physics, University of Coimbra
CITEUC
