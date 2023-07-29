# truss_builder
Calculates loads and moments on members and nodes of a user-provided truss design.

The Truss class provided by this package allows to define a set of truss nodes,
members and loads, to then generate a set of coupled equations that can be
solved by the user's solver of choice, such as numpy.linalg.lstsq. Finally, the
Truss class can use the associated results to provide the resulting loads and
moments on all members and nodes. The interface allows to specify which nodes
can slide and/or spin and which members are hinged.
