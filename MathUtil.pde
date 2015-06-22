/* Math Util 
 * provides helper methods for different math-related things
 */
 
/**
 * Bilinear Interpolation
 * Given 4 corner points, their respective values, and a point somewhere within that rectangle
 * Returns an intelligently averaged value for the middle point
 * Positions are plugged in clockwise starting with top left
 *   each position is followed by respective value
 * The desired position is plugged in last.
 * 
 * Top positions and bottom positions must have same y value respectively
 * Left positions and right positions must have same x value respectively
 */
float bilinearInterpolate(PVector topLeftPosition, float topLeftValue,
                          PVector topRightPosition, float topRightValue,
                          PVector bottomRightPosition, float bottomRightValue,
                          PVector bottomLeftPosition, float bottomLeftValue,
                          PVector desiredPosition) {
  // for sanity's sake, make sure we plugged in a rectangle
  jAssert("bilinearInterpolate: top y positions are not matching: " + topLeftPosition + " vs " + topRightPosition, topLeftPosition.y == topRightPosition.y);
  jAssert("bilinearInterpolate: bottom y positions are not matching: " + bottomLeftPosition + " vs " + bottomRightPosition, bottomLeftPosition.y == bottomRightPosition.y);
  jAssert("bilinearInterpolate: left x positions are not matching: " + topLeftPosition + " vs " + bottomLeftPosition, topLeftPosition.x == bottomLeftPosition.x);
  jAssert("bilinearInterpolate: right x positions are not matching: " + topRightPosition + " vs " + bottomRightPosition, topRightPosition.x == bottomRightPosition.x);
  
  // get the rectangle's width and height
  float rectDeltaX = topRightPosition.x - topLeftPosition.x;
  float rectDeltaY = topLeftPosition.y - bottomLeftPosition.y;
  
  // get desiredPosition's x and y fractional positions between the four corners
  float fractionX = (desiredPosition.x - topLeftPosition.x) / rectDeltaX;
  float fractionY = (desiredPosition.y - topLeftPosition.y) / rectDeltaY;
  
  if (rectDeltaX == 0)
    fractionX = 1;
  if (rectDeltaY == 0)
    fractionY = 1;
  
  // first interpolate between the top positions
  float topValue = lerp(topLeftValue, topRightValue, fractionX);
  float bottomValue = lerp(bottomLeftValue, bottomRightValue, fractionX);
  
  // now interpolate between the top and bottom
  return lerp(topValue, bottomValue, fractionY);
}

