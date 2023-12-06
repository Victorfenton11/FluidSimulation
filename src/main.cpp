#define _USE_MATH_DEFINES
#define GL_SILENCE_DEPRECATION

#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <GL/glut.h>
#include <GL/freeglut.h>
#include <cmath>
#include <unordered_map>

using namespace std;

const static glm::vec2 G = glm::vec2(0.0f, -10.0f); // Gravity vector
const static float REST_DENS = 300.0f;
const static float GAS_CONST = 2000.0f;
const static float EPS = 16.0f; // kernel radius
const static float HSQ = EPS * EPS;
const static float MASS = 2.5f; // assuming all particles have equal mass
const static float VISC = 200.0f;
const static float DT = 0.0007f;

// Gradient values
const static float POLY6 = 4.0f / (M_PI * pow(EPS, 8.0f));
const static float SPIKY_GRAD = -10.0f / (M_PI * pow(EPS, 5.0f));
const static float VISC_LAP = 40.0f / (M_PI * pow(EPS, 5.0f));

const static float BOUND_DAMPING = -0.5f;

struct Particle
{
	glm::vec2 pos, vel, force;
	float density, pressure;
	glm::vec4 color;
	Particle() : pos(0.0f), vel(0.0f), force(0.0f), density(0.0f), pressure(0.0f), color(0.0f) { }
	Particle(float x, float y) : pos(x, y), vel(0.0f), force(0.0f), density(0.0f), pressure(0.0f), color(0.0f) { }
};

struct Line
{
	glm::vec2 start, end;
	bool horizontal;

	Line() : start(0.0f), end(0.0f), horizontal(false) { }
};

static vector<Particle> particles;

const static int MAX_PARTICLES = 2500;
const static int NUM_PARTICLES = 500; // Initial number of particles for simulation: modify to view behavior with more particles

const static int WINDOW_WIDTH = 600;
const static int WINDOW_HEIGHT = 1000;
const static double VIEW_WIDTH = 1.5f * WINDOW_WIDTH;
const static double VIEW_HEIGHT = 1.5f * WINDOW_HEIGHT;

// Spatial Hashing: Note that this increases the speed of the simulation significantly
static unordered_map<int, vector<Particle>> SHGrid;
int cellDim = VIEW_WIDTH / 60;

// Maze params
int mazeLength;
float mazeThreshold;
static vector<Line> mazeWalls;
const static float marginTop = VIEW_HEIGHT - 400.0f;
static bool particleOnly = false;

// Function prototypes
void initGLUT(int*, char**);
void display(void);
void reshape(GLint, GLint);
void keyPress(unsigned char, int, int);
void keyRelease(unsigned char, int, int);
void updateParticles(void);
void spatialHashing(void);
int SHhash(Particle);
vector<Particle> findNeighbors(Particle);
void computeDensityPressure(void);
void computeExtForces(void);
void moveParticles(void);
void initSPHSystem(void);
void drawMaze();


int main(int argc, char** argv) {
	cout << "CS334 Final Project: Falling Water - Victor Fenton Aguilar" << endl << endl;

	try {
		// Create the window and menu
		initGLUT(&argc, argv);
	}
	catch (const exception& e) {
		// Handle any errors
		cerr << "Fatal error: " << e.what() << endl;
		return -1;
	}

	do
	{
		cout << "Enter desired maze dimension (minimum 5, maximum 50): ";
		scanf("%d", &mazeLength);
		cout << "Enter a probability threshold for the random generation of maze walls (minimum 0.3, maximum 0.7): ";
		scanf("%f", &mazeThreshold);
	} while (mazeLength < 5 || mazeLength > 50 || mazeThreshold < 0.3 || mazeThreshold > 0.7);

	cout << endl;
	drawMaze();

	initSPHSystem();
	cout << endl << "Keyboard Controls:" << endl;
	cout << "\tESC to exit." << endl;
	cout << "\tSpace to add more particles to the simulation." << endl;
	cout << "\t\'r\' to restart the flow of particles." << endl;
	cout << "\t\'m\' to generate a new random maze and restart the SPH simulation." << endl;
	cout << "\t\'p\' to restart the simulation without the maze, just fluid particles." << endl << endl;

	glutMainLoop();

	return 0;
}


void initGLUT(int* argc, char** argv) {
	glutInit(argc, argv);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
	glutInitWindowPosition(WINDOW_WIDTH, 0);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB);
	glutCreateWindow("Particle Fluid Simulation | Victor Fenton Aguilar");

	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardUpFunc(keyRelease);
	glutKeyboardFunc(keyPress);
	glutIdleFunc(updateParticles);

	glClearColor(0.0f, 0.0f, 0.0f, 0.5f);
	glEnable(GL_POINT_SMOOTH);
	glPointSize(EPS / 2.0f);
	glMatrixMode(GL_PROJECTION);
}


void display(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glLoadIdentity();
	glOrtho(0, VIEW_WIDTH, 0, VIEW_HEIGHT, 0, 1);

	glBegin(GL_POINTS);
	for (auto& p : particles)
	{
		glColor4f(p.color.x, p.color.y, p.color.z, p.color.w);
		glVertex2f(p.pos.x, p.pos.y);
	}
	glEnd();

	glBegin(GL_LINES);
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	glLineWidth(5.0f);
	for (auto& l : mazeWalls)
	{
		glVertex2f(l.start.x, l.start.y);
		glVertex2f(l.end.x, l.end.y);
	}
	glEnd();

	glutSwapBuffers();
}

void drawMaze()
{
	float segmentLength;
	
	segmentLength = VIEW_WIDTH / mazeLength;

	Line left = Line();
	left.horizontal = true;
	left.start = glm::vec2(0.0f, marginTop);
	float lStart = ((mazeLength / 2) - (mazeLength % 2 == 0))* segmentLength;
	left.end = glm::vec2(lStart, marginTop); 
	mazeWalls.push_back(left);
	Line right = Line();
	right.horizontal = true;
	float rStart = ((mazeLength / 2) + 1) * segmentLength;
	right.start = glm::vec2(rStart, marginTop);
	right.end = glm::vec2(VIEW_WIDTH, marginTop);
	mazeWalls.push_back(right);

	float random;
	for (int i = 0; i < mazeLength; i++)
	{
		for (int j = 1; j <= mazeLength; j++)
		{
			if (j != mazeLength)
			{
				random = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
				if (random > mazeThreshold + 0.1f) {
					Line lh = Line();
					lh.start = glm::vec2(segmentLength * i, marginTop - (segmentLength * j));
					lh.end = lh.start + glm::vec2(segmentLength, 0.0f);
					lh.horizontal = true;
					mazeWalls.push_back(lh);
				}
			}

			if (i != mazeLength - 1)
			{
				random = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
				if (random > mazeThreshold)
				{
					Line lv = Line();
					lv.start = glm::vec2(segmentLength * (i + 1), marginTop - (segmentLength * j));
					lv.end = lv.start + glm::vec2(0.0f, segmentLength);
					lv.horizontal = false;
					mazeWalls.push_back(lv);
				}
			}
		}
	}

	Line bottomLeft = Line();
	bottomLeft.horizontal = true;
	float bh = marginTop - (mazeLength * segmentLength);
	bottomLeft.start = glm::vec2(0.0f, bh);
	bottomLeft.end = glm::vec2(lStart, bh);
	mazeWalls.push_back(bottomLeft);
	Line bottomRight = Line();
	bottomRight.horizontal = true;
	bottomRight.start = glm::vec2(rStart, bh);
	bottomRight.end = glm::vec2(VIEW_WIDTH, bh);
	mazeWalls.push_back(bottomRight);
}


void initSPHSystem(void)
{
	cout << "Starting \'dam break\' with " << NUM_PARTICLES << " particles" << endl;
	for (float y = marginTop + 4*EPS; y < VIEW_HEIGHT - EPS; y += EPS)
	{
		for (float x = VIEW_WIDTH / 3; x <= 2*VIEW_WIDTH / 3; x += EPS)
		{
			if (particles.size() < NUM_PARTICLES)
			{
				Particle p = Particle();
				float offset = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
				if (offset < 0.1) p.color = glm::vec4(0.7f, 0.9f, 0.98f, 1.0f);
				else if (offset < 0.7) p.color = glm::vec4(0.2f, 0.365f, 0.655f, 1.0f);
				else p.color = glm::vec4(0.23f, 0.54f, 0.77f, 1.0f);
				p.pos = glm::vec2(x + offset, y);
				particles.push_back(p);
			}
			else
			{
				return;
			}
		}
	}
}


void updateParticles(void)
{
	spatialHashing();
	computeDensityPressure();
	computeExtForces();
	moveParticles();

	glutPostRedisplay();
}

void spatialHashing()
{
	SHGrid.clear();
	for (auto &p : particles)
	{
		int cell = SHhash(p);
		SHGrid[cell].push_back(p);
	}
}


int SHhash(Particle p)
{
	return static_cast<int>(p.pos.x / cellDim) + (static_cast<int>(p.pos.y / cellDim) * 1000);
}

vector<Particle> findNeighbors(Particle p)
{
	vector<Particle> neighbors;

	int cell = SHhash(p);
	for (int i = -1; i <= 1; i++)
	{
		for (int j = -1; j <= 1; j++)
		{
			int key = cell + (i * 1000) + j;
			if (SHGrid.find(key) == SHGrid.end()) continue;
			for (auto& pn : SHGrid[key])
			{
				neighbors.push_back(pn);
			}
		}
	}
	return neighbors;
}

void computeDensityPressure(void)
{
	for (auto& p : particles)
	{
		p.density = 0.f;
		for (auto& pn : findNeighbors(p))
		{
			glm::vec2 r = pn.pos - p.pos;
			float r2 = glm::pow(glm::length(r), 2.0);

			if (r2 < HSQ)
			{
				p.density += MASS * POLY6 * pow(HSQ - r2, 3.f);
			}
		}
		p.pressure = GAS_CONST * (p.density - REST_DENS);
	}
}


void computeExtForces(void)
{
	for (auto& pi : particles)
	{
		glm::vec2 fpress(0.f, 0.f);
		glm::vec2 fvisc(0.f, 0.f);
		for (auto &pj : particles)
		{
			if (&pi == &pj) continue;
			glm::vec2 rij = pj.pos - pi.pos;
			float r = glm::length(rij);

			if (r < EPS)
			{
				fpress += -glm::normalize(rij) * MASS * (pi.pressure + pj.pressure) / (2.f * pj.density) * SPIKY_GRAD * pow(EPS - r, 3.f);
				fvisc += VISC * MASS * (pj.vel - pi.vel) / pj.density * VISC_LAP * (EPS - r);
			}
		}
		glm::vec2 fgrav = G * MASS / pi.density;
		pi.force = fpress + fvisc + fgrav;
	}
}


void moveParticles(void)
{
	for (auto& p : particles)
	{
		p.vel += DT * p.force / p.density;
		if (glm::length(p.vel) < glm::pow(10, -5)) p.vel = glm::vec2(0.0f);
		p.pos += DT * p.vel;

		if (p.pos.x - EPS < 0.f)
		{
			p.vel.x *= BOUND_DAMPING;
			p.pos.x = EPS;
		}
		if (p.pos.x + EPS > VIEW_WIDTH)
		{
			p.vel.x *= BOUND_DAMPING;
			p.pos.x = VIEW_WIDTH - EPS;
		}
		if (p.pos.y - EPS < 0.f)
		{
			p.vel.y *= BOUND_DAMPING;
			p.pos.y = EPS;
			if (!particleOnly)
			{
				cout << endl << "This Random Maze is Solvable!!" << endl;
				particleOnly = true;
			}
		}
		if (p.pos.y + EPS > VIEW_HEIGHT)
		{
			p.vel.y *= BOUND_DAMPING;
			p.pos.y = VIEW_HEIGHT - EPS;
		}
		for (auto& line : mazeWalls)
		{
			float radius = EPS / 2.0f;
			if (!line.horizontal)
			{
				if ((glm::abs(p.pos.x - line.start.x) <= radius) && (line.start.y <= p.pos.y && p.pos.y <= line.end.y))
				{
					p.vel.x *= BOUND_DAMPING;
					p.pos.x += p.pos.x - line.start.x;
				}
			}
			else {
				if ((glm::abs(p.pos.y - line.start.y) <= radius) && (line.start.x <= p.pos.x && p.pos.x <= line.end.x))
				{
					p.vel.y *= BOUND_DAMPING;
					p.pos.y += p.pos.y - line.start.y;
				}
			}
		}
	}
}


void keyPress(unsigned char key, int x, int y)
{
	switch (key)
	{
	case ' ':
		if (particles.size() >= MAX_PARTICLES)
		{
			std::cout << "Maximum allowed number of particles reached" << std::endl;
		}
		else
		{
			for (float y = marginTop + 4*EPS; y < VIEW_HEIGHT - EPS; y += EPS * 0.95f)
			{
				for (float x = VIEW_WIDTH / 3; x <= 2 * VIEW_WIDTH / 3; x += EPS)
				{
					if (particles.size() < MAX_PARTICLES)
					{
						Particle p = Particle();
						float offset = static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
						if (offset < 0.1) p.color = glm::vec4(0.7f, 0.9f, 0.98f, 1.0f);
						else if (offset < 0.7) p.color = glm::vec4(0.2f, 0.365f, 0.655f, 1.0f);
						else p.color = glm::vec4(0.23f, 0.54f, 0.77f, 1.0f);
						p.pos = glm::vec2(x + offset, y);
						particles.push_back(p);
					}
				}
			}
		}
		break;
	case 'r':
	case 'R':
		particles.clear();
		initSPHSystem();
		break;
	case 'm':
		particles.clear();
		mazeWalls.clear();
		drawMaze();
		initSPHSystem();
		particleOnly = false;
		break;
	case 'p':
		particles.clear();
		mazeWalls.clear();
		initSPHSystem();
		particleOnly = true;
		break;
	}
}

// Called when a key is released
void keyRelease(unsigned char key, int x, int y) {
	switch (key) {
	case 27:	// Escape key
		glutLeaveMainLoop();	
		break;
	}
}

// Called when the window is resized
void reshape(GLint width, GLint height) {
	glViewport(0, 0, width, height);
	return;
}