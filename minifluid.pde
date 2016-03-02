FluidGrid fluidGrid;

int cellWidth = 35;
void setup() {
  // 1280 720
  size(1280, 720, P2D);
  
  fluidGrid = new FluidGrid(1, 2, 2);
//  noLoop();
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
