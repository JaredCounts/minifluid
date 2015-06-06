/**
 * Particle
 * Lives in fluid region and follows velocity fields
 */
class Particle {
  // position and velocity vectors
  // only to be modified by particle
  private final PVector position, velocity;
  
  Particle(PVector position, PVector velocity) {
    // to prevent rep exposure, we copy the given position and velocities
    this.position = position.get();
    this.velocity = velocity.get();
  }
}
