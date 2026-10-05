# Lateral flows

!!! warning "In development"
    Runoff routing is in an open pull request,
    [#33](https://github.com/BiophysicalEcology/MicroclimateMapper.jl/pull/33), and not yet in a release. This page
    describes its design. Cold-air drainage is planned.

Other than the horizon angles, Microclimate.jl, is not "pixel-aware". It treats each point as part of an infinite, 
uniform plane: nothing flows in from the side. In reality, water and cold air may run downhill and pool in the 
lowlands. Those flows need the grid context and are thus the domain of MicroclimateMapper.

## Runoff routing

Each cell already computes the surface water that its soil cannot take, `runoff_generated`. Routing passes it
downhill as inflow to the next cell's water balance.

```julia
model = MicroMapModel(; micro_model, dem_source = SRTM, weather_source = NCEP{SurfaceFlux, 1},
    routing_model = SurfaceRunoffRouting(; sinks = TerminalSinks()))
```

Routing works on grids, not points. It adds a `runoff_generated` layer to the output.

The design keeps the speed of independent cells:

- **A flow graph built once.** From the run grid's DEM, each cell drains to its steepest downhill neighbour of
  eight (D8). The cells form a forest of trees whose roots are outlets at the edge of the grid or closed
  depressions.
- **Cells in dependency order, not time steps.** A cell's inflow is the whole time series of its upstream cells'
  outflow, which is then a forcing like the weather. So each cell still solves its full period in one call, as
  without routing. There is no global time loop.
- **In parallel.** A cell is ready once all its upstream cells are done. Threads take ready cells from a shared
  queue and follow each flow path downstream. Inflows are summed in a fixed order, so the results do not depend on
  the number of threads.
- **Depressions** are either terminal sinks, where water collects and leaves only by evaporation and infiltration
  (`TerminalSinks()`, the default), or assumed full and spilling at their sill (`SpillOver()`).

The added cost is small, because each cell is still solved once.

## Cold-air drainage

On clear, still nights, air cooled at the ground flows downslope and pools in hollows and valleys, which can be
several degrees colder than the slopes above. The same flow graph is a candidate for routing it.

## Limitations

This approach is fast but it results in the final basin pixel taking all the residual flow to form a big column of 
air or water. Ideally the pooling fluid would fill the basin but such a calculation would need a different, slower
algorithm and is not presently planned for this package.

## Beyond routing

Routing along a fixed graph suits flows that go one way, downhill. Flows that need a lateral solve, such as
groundwater, or air moving in both directions, need a model built on a spatial discretisation.
[ClimaLand.jl](https://clima.github.io/ClimaLand.jl/stable/), the land model of the CliMA project, builds on
ClimaCore, whose spaces handle lateral flows. Coupling to it is under consideration. Microclimate.jl's
[Introduction](https://biophysicalecology.github.io/Microclimate.jl/dev/manual/introduction) compares it with
the point model.
