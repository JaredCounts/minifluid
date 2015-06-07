/**
 * Particle
 * Lives in fluid region and follows velocity fields
 * Particles are massless, so they have no physical effect on the simulation
 * They are simply used for visualizing the liquid!
 */
 
class Particle {
  // position and velocity vectors
  // only to be modified by particle
  private final PVector position;
  
  // a reference to the shared fluid grid among particles
  private final FluidGrid fluidGrid;
  
  /**
   * Construct a particle
   * position must be within the fluid grid's region
   */
  Particle(FluidGrid fluidGrid, PVector position) {
    assert(position.x < fluidGrid.getRegionWidth());
    assert(position.y < fluidGrid.getRegionHeight());
    
    this.fluidGrid = fluidGrid;
    
    // to prevent rep exposure, we copy the given position
    this.position = position.get();
  }
  
  /**
   * update(elapsedTime)
   * Integrate's the particle's position through fluidGrid's velocity field
   * @param elapsedTime Time elapsed since the last update.
   */
  void update(float elapsedTime) {
    PVector velocity = fluidGrid.getVelocityAt(position);
    velocity.mult(elapsedTime);
    position.add(velocity);
  }
}
