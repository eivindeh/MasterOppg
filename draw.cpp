#include <stdio.h> // printf, sprintf, fflush, stdout
#include <math.h>  // atan, pow, ceil (for std_easy_font.h)

// OpenGL/GLFW stuff
#include <GLFW/glfw3.h>

GLFWwindow* window;
double resx = 1440/2, resy = 1150/2;
double prevx = -1, prevy = -1;     // to track how much the mouse moved between frames
double cam_x = -0.1, cam_y = -0.1; // world coordinates of lower-left corner of the window
double cam_height = 1.2;
double cam_width = cam_height*resx/double(resy); 
int clickedButtons = 0;

enum buttonMaps { FIRST_BUTTON=1, SECOND_BUTTON=2, THIRD_BUTTON=4, FOURTH_BUTTON=8, FIFTH_BUTTON=16, NO_BUTTON=0 };
enum modifierMaps { CTRL=2, SHIFT=1, ALT=4, META=8, NO_MODIFIER=0 };

void initGL(void);
void windowsize_callback(GLFWwindow*, int, int);
void key_callback(GLFWwindow*, int, int, int, int);
void mousebutton_callback(GLFWwindow*, int, int, int);
void mousepos_callback(GLFWwindow*, double, double);
void mousewheel_callback(GLFWwindow*, double, double);
void Draw(void);


// Text stuff
#include "stb_easy_font.h"

enum TextAlignment {ALIGN_MID, ALIGN_LEFT, ALIGN_RIGHT};

void print_string(float, float, char*, float, float, float);
void print_string_world(float, float, char, float, float, float, float, int);
void print_string_screen(float, float, char*, float, float, float, float);

bool PAUSE = true;

int mainLoop() {
	cam_height = yLen*2;
	cam_width = cam_height*resx/resy;
	initGL();
	printf("%s\n", glGetString(GL_VERSION)); fflush(stdout);
	// Main loop
	while ( !glfwWindowShouldClose(window)) {	
		if (!PAUSE) {
			int iter = 1;
			if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS) {
				iter *= 10;
			}
			if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS) {
				iter *= 10;
			}
			for (int i = 0; i < iter; i++)
					timestep();
		}
		Draw();
		glfwPollEvents();
	}

	return 0;
}



// Printing text using an 8 pixel bitmap font (amiga clone?)
// Assume ortogonal projection and screen coordinates
// Draws at depth 0
// from https://github.com/nothings/stb/
void print_string(float x, float y, char *text, float r, float g, float b) {
	static char buffer[99999]; // ~500 chars, enough to fill screen
	int num_quads;
	num_quads = stb_easy_font_print(x, y, text, NULL, buffer, sizeof(buffer));
	glColor3f(r,g,b);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(2, GL_FLOAT, 16, buffer);
	glDrawArrays(GL_QUADS, 0, num_quads*4);
	glDisableClientState(GL_VERTEX_ARRAY);
}

// Draw text in screen coordinates
// Saves matrix states
void print_string_screen(float x, float y, char *text, float r, float g, float b, float scale = 1.0) {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0, resx, resy, 0.0, -1.0, 1.0); // top to bottom for text

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();
	glTranslatef(0.0, 0.0, 0.1);    // Move text layer up 0.1 units to be above 0-level
	glScalef(scale, scale, scale);

	print_string(x, y, text, r, g, b);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();
}

// Draw text in world coordinates
// Saves matrix states
// Uses an arbitrary scale
// ALIGN_MID (default), ALIGN_LEFT or ALIGN_RIGHT to justify/align the text relative to position
void print_string_world(float x, float y, char *text, float r, float g, float b, float scale, int alignment=ALIGN_MID) {
	glMatrixMode(GL_PROJECTION);
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0, resx, resy, 0.0, -1000, 1000);

	glMatrixMode(GL_MODELVIEW);
	glPushMatrix();

	int text_width = stb_easy_font_width(text);
	scale *= 0.01/(cam_width/resy); // convert scale from world to screen coordinates using camera information
	
	if (alignment == ALIGN_MID) {
		glTranslatef((x - cam_x)/cam_width * resx - 0.5*(text_width-1.0)*scale, (1-(y - cam_y)/cam_height)*resy - (8-1)*scale/2 , 0);
	} else if (alignment == ALIGN_RIGHT) {
		glTranslatef((x - cam_x)/cam_width * resx - 1.0*(text_width-1.0)*scale, (1-(y - cam_y)/cam_height)*resy - (8-1)*scale/2 , 0);
	} else if (alignment == ALIGN_LEFT) {
		glTranslatef((x - cam_x)/cam_width * resx - 0.0*(text_width-1.0)*scale, (1-(y - cam_y)/cam_height)*resy - (8-1)*scale/2 , 0);
	}
	glScalef(scale, scale, scale);
	print_string(0.0, 0.0, text, r, g, b);

	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	glMatrixMode(GL_MODELVIEW);
	glPopMatrix();

}

struct vec3 {
	double x, y, z;

	vec3(double x, double y, double z) : x(x), y(y), z(z) {}
};

void drawBubble(vec3 color, double x1, double x2, double y1, double y2, double s, double t, double len) {
    double theta = acos((x2 - x1)/sqrt((x2 - x1)*(x2 - x1) + (y2 - y1)*(y2 - y1)));

    double X1 = x1 + len/2*(x2 - x1) + (x2-x1)*(1.0 - len)*s;
    double X2 = x2 - len/2*(x2 - x1) - (x2-x1)*(1.0 - len)*(1-t);

    double Y1 = y1 + len/2*(y2 - y1) + (y2-y1)*(1.0 - len)*s;
    double Y2 = y2 - len/2*(y2 - y1) - (y2-y1)*(1.0 - len)*(1-t);

    double L = sqrt(((y2 - len/2*(y2 - y1)) - (y1 + len/2*(y2 - y1)))*((y2 - len/2*(y2 - y1))- (y1 + len/2*(y2 - y1))) + ((x2 - len/2*(x2 - x1)) - (x1 + len/2*(x2 - x1)))*((x2 - len/2*(x2 - x1)) - (x1 + len/2*(x2 - x1))));
    len = len*L;

    glColor3f(color.x, color.y, color.z);
    glBegin(GL_QUADS);
    glVertex2f(X1 - len*sin(theta), Y1 + len*cos(theta));
    glVertex2f(X1 + len*sin(theta), Y1 - len*cos(theta));
    glVertex2f(X2 + len*sin(theta), Y2 - len*cos(theta));
    glVertex2f(X2 - len*sin(theta), Y2 + len*cos(theta));
    glEnd();
}

// Main drawing function
void Draw() {

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Draw non-text stuff
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glOrtho(cam_x, cam_x + cam_width, cam_y, cam_y + cam_height, -1.0, 1.0);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();


	// Draw coordinate system lines
	glColor3f(0.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(cam_x,             0.0,               0.05);
	glVertex3f(cam_x + cam_width, 0.0,               0.05);
	glVertex3f(0.0,               cam_y,              0.05);
	glVertex3f(0.0,               cam_y + cam_height, 0.05);
	glEnd();

	// Draw coordinate numbers
	static char str_coordinate[12];
	for (int i = cam_x - 1.0; i < cam_x + cam_width + 1.0; i++) {
		if (i == 0)
			continue;
		sprintf(str_coordinate, "%d", i);
		print_string_world(i, -0.05, str_coordinate, 0, 0, 0, 1.0);
	}
	for (int j = cam_y - 1.0; j < cam_y + cam_height + 1.0; j++) {
		sprintf(str_coordinate, "%d", j);
		print_string_world(-0.01, (j > 0 ? j : j - 0.05), str_coordinate, 0.0, 0.0, 0.0, 1.0, ALIGN_RIGHT);
	}

	double dy = (sin(PI/6) + 1)/cos(PI/6);

	 for (int i = 1; i < nPipes;  i++) {
                double x1 = pipes[i].x1;
                double y1 = pipes[i].y1;
                double x2 = pipes[i].x2;
                double y2 = pipes[i].y2;
                

                if (x2 - x1 > xLen/2) {
                	x1 += xLen;
                } else if (x2 - x1 < -xLen/2) {
                	x2 += xLen;
                }

                if (y2 - y1 < -yLen/2) {
                	y2 += yLen*dy;
                }


                double len = linkRadius[l];

                drawBubble(vec3(0,0,1), x1, x2, y1, y2, 0.0, 1.0, len);
                glTranslatef(0.0, 0.0, 0.1);
                for (int b = 0; b < linkNumBubbles[l]; b++) {
                    drawBubble(vec3(1,0,0), x1, x2, y1, y2, bubbleStart[l][b], bubbleStop[l][b], len);
                }
                glTranslatef(0.0, 0.0, -0.1);
            }
        }
    }



	// Draw Text
	static char str_text[256];
	static int frame = 0;
	int pos = 0;
	sprintf(str_text, "This is frame #%d", frame++);
	print_string_screen(3, 3 + 10*pos++, str_text, 0.0, 0.0, 0.0, 1.5);
	sprintf(str_text, "This is on the next line. Camera position = (%.4f, %.4f), (%.4f, %.4f)", cam_x, cam_y, cam_x + cam_width, cam_y + cam_height);
	print_string_screen(3, 3 + 10*pos++, str_text, 0.0, 0.0, 0.0, 1.5);

	// done, possibly wait for vsync
	glfwSwapBuffers(window);
}

// Initialize GLFW and OpenGL:
// Craetes a window, make it the current content
// Setup orthogonal coordinates
// Set background to white
// Tries to enable VSync
// Set up callback functions
void initGL() {
	printf("Initializing OpenGL/GLFW\n"); 
	if (!glfwInit()) {
		printf("Could not initialize\n");
		exit(-1);
	}

	glfwWindowHint(GLFW_SAMPLES , 4);
	window = glfwCreateWindow(resx, resy, "window title", 0, 0);
	if (!window) {
		printf("Could not open glfw window\n");
		glfwTerminate();
		exit(-2);
	}
	glfwMakeContextCurrent(window);

	glLoadIdentity();
	glOrtho(cam_x, cam_x + cam_width, cam_y, cam_y + cam_height, -1.0, 1.0);	
	glClearColor(1.0, 1.0, 1.0, 1.0);
	glfwSwapInterval(1);
	glEnable(GL_DEPTH_TEST);

	glfwSetKeyCallback(window, key_callback);
	glfwSetMouseButtonCallback(window, mousebutton_callback);
	glfwSetScrollCallback(window, mousewheel_callback);
	glfwSetCursorPosCallback(window, mousepos_callback);
	glfwSetWindowSizeCallback(window, windowsize_callback);
}

// Callback function called every time the window size change
// Adjusts the camera width and heigh so that the scale stays the same
// Resets projection matrix
void windowsize_callback(GLFWwindow* /*win*/, int width, int height) {
	double distance_per_pixel = cam_height/resy; // assuming cam_height/resy == cam_width/resx

	resx = width;
	resy = height;
	cam_width = distance_per_pixel*resx;
	cam_height = distance_per_pixel*resy;

	glLoadIdentity();
	glViewport(0, 0, resx, resy);
	glOrtho(cam_x, cam_x + cam_width, cam_y, cam_y + cam_height, -1, 1);
}

// Callback function called every time a keyboard key is pressed, released or held down
void key_callback(GLFWwindow* win, int key, int /*scancode*/, int action, int /*mods*/) {
	//printf("key = %d, scancode = %d, action = %d, mods = %d\n", key, scancode, action, mods); fflush(stdout);

	// Close window if escape is released
	if (key == GLFW_KEY_ESCAPE && action == GLFW_RELEASE) {
		glfwSetWindowShouldClose(win, GL_TRUE);
	}

	if (key == GLFW_KEY_SPACE && action) {
		timestep();
	}

	if (key == GLFW_KEY_P && action) {
		PAUSE = !PAUSE;
	}

}

// Callback function called every time a mouse button pressed or released
void mousebutton_callback(GLFWwindow* win, int button, int action, int /*mods*/) {
	// get current cursor position, convert to world coordinates
	glfwGetCursorPos(win, &prevx, &prevy);
	//double xend = cam_x + prevx*(cam_width)/resx;
	//double yend = cam_y + (resy - prevy)*cam_height/resy;	

	//printf("button = %d, action = %d, mods = %d at (%f %f)\n", button, action, mods, xend, yend); fflush(stdout);

	// To track the state of buttons outside this function
	if (action == 1)
		clickedButtons |= (1 << button);
	else
		clickedButtons &= ~(1 << button);


	// Test each button
	if (clickedButtons&FIRST_BUTTON) {
		
	} else if (clickedButtons&SECOND_BUTTON) {

	} else if (clickedButtons&THIRD_BUTTON) {

	} else if (clickedButtons&FOURTH_BUTTON) {

	} else if (clickedButtons&FIFTH_BUTTON) {

	}
}

// Callback function called every time a the mouse is moved
void mousepos_callback(GLFWwindow* /*win*/, double xpos, double ypos) {
	// move the camera if the first mouse button is held down
	// the cursor will stay on at the same location relative to world coordinates after movement
	if (clickedButtons&FIRST_BUTTON) {
		cam_x -= (xpos-prevx)/resx*cam_width;
		cam_y += (ypos-prevy)/resy*cam_height;

		prevx = xpos;
		prevy = ypos;
	} else if (clickedButtons&SECOND_BUTTON) {

	} else if (clickedButtons&THIRD_BUTTON) {

	} else if (clickedButtons&FOURTH_BUTTON) {

	} else if (clickedButtons&FIFTH_BUTTON) {

	}
}

// Callback function called every time a the mouse scroll wheel is moved
// yoffset = up-down
// xoffset = left-right
void mousewheel_callback(GLFWwindow* win, double /*xoffset*/, double yoffset) {
	double zoomFactor = pow(0.95,yoffset);

	glfwGetCursorPos(win, &prevx, &prevy);                  // get mouse position in window
	double xend = cam_x + prevx*(cam_width)/resx;			// convert screen position to world cooridnates
	double yend = cam_y + (resy - prevy)*cam_height/resy;	
	cam_x = (1.0-zoomFactor)*xend + zoomFactor*cam_x;		// update lower left corner
	cam_y = (1.0-zoomFactor)*yend + zoomFactor*cam_y;

	cam_width *= zoomFactor;
	cam_height *= zoomFactor;

	// zoom the camera by keeping the center constant and move the edges accordingly
	//cam_x = cam_x + cam_width*(1.0 - zoomFactor)/2.0;
	//cam_y = cam_y + cam_height*(1.0 - zoomFactor)/2.0;

	// cam_width *= zoomFactor;
	//cam_height *= zoomFactor;
}
