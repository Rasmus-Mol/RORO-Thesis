# RoRo Vessel Optimizer

## Overview
The RoRo Vessel Optimizer is a Julia-based operations research project that utilizes the JuMP package to model and optimize the operations of Roll-on/Roll-off (RoRo) vessels. The project aims to enhance the efficiency of cargo loading and unloading processes while adhering to operational constraints.

## Project Structure
- `src/StowagePlanner.jl`: Main module for the application.
- `src/model/base_model.jl`: Defines the base model structure, including the objective function and decision variables.
- `src/model/stability.jl`: Contains stability constraints for the vessel.
- `src/model/hazardous.jl`: Contains hazardous cargo handling constraints.
- `src/representation/`: Contains data structures for vessel, cargo, and slots.
- `test/runtests.jl`: Runs all tests in the project.

## Setup Instructions
1. Ensure you have Julia installed on your machine. If not, please use juliaup and julia v1.11.3
2. Clone the repository:
3. Navigate to the project directory:
4. Install the required dependencies:
```julia
using Pkg
Pkg.instantiate()

## Usage
To run the optimization model, execute the following command in the Julia REPL:
```bash
julia main.jl
```
