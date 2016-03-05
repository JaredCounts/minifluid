FluidGrid fluidGrid;

int cellWidth = 35;
void setup() {
  // 1280 720
  size(1280, 720, P2D);
  
  fluidGrid = new FluidGrid(cellWidth, width, height);
//  noLoop();
  //frameRate(2);
}
void draw() {
  // physics
  fluidGrid.solve();
  
  background(255);
  fluidGrid.draw();
}

void mouseDragged() {
  int cellX = mouseX / cellWidth;
  int cellY = mouseY / cellWidth;
  
  float mouseSpeedX = mouseX - pmouseX;
  float mouseSpeedY = mouseY - pmouseY;
  //println(mouseSpeedX);
  
  fluidGrid.cells[cellX][cellY].velocityXLeft += mouseSpeedX;
  fluidGrid.cells[cellX][cellY].velocityYTop += mouseSpeedY;
  fluidGrid.cells[cellX][cellY].pressure += 1;
  
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