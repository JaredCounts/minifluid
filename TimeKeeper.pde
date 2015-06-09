/**
 * Time Keeper
 * Mutable singleton class
 * Processing doesn't like static classes, so I need to just make one here
 * For keeping track of how much time we need to solve for per each frame
 */
TimeKeeper timeKeeper = new TimeKeeper();
class TimeKeeper {
  private long lastSolveTime;
  private long leftOverSolveTime;
  private boolean firstSolve = true;
  
  /**
   * Get amount of time elapsed since last call in seconds
   * plus amount of "left over" time that couldn't be solved for during the last frame
   */
  public float getTimeToSolveFor() {
    float secondsToSolveFor;
    if (!firstSolve) {
      // compute amount of time to solve for
      long currentMillis = millis();
      long elapsedMillis = currentMillis - lastSolveTime;
      long totalMillis = elapsedMillis + leftOverSolveTime;
      secondsToSolveFor = totalMillis / 1000.f;
    }
    else { // first solve
      secondsToSolveFor = 1.f / 60.f; // solve as if we've been running at 60 fps
      firstSolve = false;
    }
    
    // reset time keeper variables
    lastSolveTime = millis();
    leftOverSolveTime = 0;
    
    return secondsToSolveFor;
  }
  
  /**
   * Set amount of time to be added to the time to solve for for the next frame.
   */
  public void setLeftOverSolveTime(float seconds) {
    long millis = (long)(seconds * 1000);
    leftOverSolveTime = millis;
  }
}
