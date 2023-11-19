# Computed Fluid Dynamic Rust Wasm lib

##Description
This is a rust CFD lib created with the intent to be used in a web enviroment. This Computed Fluid Dynamic simulation uses Euler based fluid method. 
The model makes the following assumptions: 
1. Incompressible flow
2. Inviscid Fluid
3. The motion of a particle is entirely 2-Dimentional

NOTE: this code does not take care of applying boundary conditions. The purpose for this code is to give a base to work with the following conditions:
1. Wall conditions (obstacles)
2. Heat transfer(materials heat transfer characteristics)
3. etc ... (any applicable condition to your case)
 
has to be implemented in the respective simulation rendering env.
