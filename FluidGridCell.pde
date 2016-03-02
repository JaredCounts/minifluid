/** 
 * Fluid Grid Cell
 * Mutable object which acts as a container for pressure, and has references to its edge velocities
 */
class FluidGridCell {
  // pressure at the center of this cell
  float pressure; 

  // velocity objects which are shared with adjacent cells
  // must set these manually
  Float velocityXLeft, velocityXRight, velocityYTop, velocityYBottom;

  // must set positions manually
  // position is the center of the cell
  PVector position;
  // each edge position is the center of each edge of the cell
  PVector topEdgePosition, rightEdgePosition, bottomEdgePosition, leftEdgePosition;
  
  // flag for whether this cell is filled with liquid
  boolean hasLiquid;

  /**
   * Construct an empty FluidGridCell
   * note that edge velocities must be set manually. 
   */
  FluidGridCell() {
    pressure = random(100);
    hasLiquid = false;
  }
  
  void draw() {
    PVector velocity = new PVector((velocityXLeft + velocityXRight) / 2.f,
                                   (velocityYTop + velocityYBottom) / 2.f);
    
    stroke(0);
    line(position.x, position.y, // from
         position.x + velocity.x, position.y + velocity.y); // to
  }
  
  PVector getCenterVelocity() {
    return new PVector((velocityXLeft + velocityXRight) / 2.f,
                                   (velocityYTop + velocityYBottom) / 2.f);
  }
  
}
