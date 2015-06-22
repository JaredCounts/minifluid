FluidGrid fluidGrid;

void setup() {
  size(1280, 720, P2D);
  
  fluidGrid = new FluidGrid(10, 1280, 720);
}
void draw() {
  // physics
  fluidGrid.solve();
  
  background(255);
  fluidGrid.draw();
}

// our own assert because Processing's assertions suck.
void jAssert(String message, boolean assertion) {
  if (!assertion) {
    println();
    println("Assertion error");
    StackTraceElement[] stackTrace = Thread.currentThread().getStackTrace();
    for (int i = 2; i < stackTrace.length; i++)
      println("\t" + stackTrace[i]);
    println("Error message");
    println("\t" + message);
    println();
    assert(assertion);
  } 
}
