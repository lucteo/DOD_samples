#include <math.h>
#include <vector>

#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

using namespace std;

static const int g_numPoints = 10 * 1000 * 1000;

static const float g_worldSize = 300 * 1000.0f;
static const int g_viewSize = 800;
static const float g_maxVelocity = 3.0f;

namespace Util
{
    ////////////////////////////////////////////////////////////////////////////
    struct Vec2 {
        float x;
        float y;

        Vec2() : x(0.0f), y(0.0f) {}
        Vec2(float xx, float yy) : x(xx), y(yy) {}
    };

    float randBetween(float minVal, float maxVal)
    {
        return minVal + (maxVal - minVal) * (static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
    }

    struct RotMatrix {
        float cos_theta;
        float sin_theta;
    };

    RotMatrix getRotMatrix(float angle)
    {
        return RotMatrix{ cos(angle), sin(angle) };
    }

    Vec2 rotate(RotMatrix mat, Vec2 vec)
    {
        return Vec2(vec.x * mat.cos_theta - vec.y * mat.sin_theta,
                    vec.x * mat.sin_theta + vec.y * mat.cos_theta);
    }

    void drawPoint(Vec2 pt)
    {
        if ( pt.x < g_viewSize && pt.y < g_viewSize )
            glVertex2f(pt.x * 2.0f / g_viewSize - 1.0f, pt.y * 2.0f / g_viewSize - 1.0f);
    }

    ////////////////////////////////////////////////////////////////////////////
    class Model {
    };
}
using namespace Util;

namespace ClassicOOP
{
    ////////////////////////////////////////////////////////////////////////////
    class IGameObject {
    public:
        virtual ~IGameObject() {}
        virtual void advance(RotMatrix velocityRotation) = 0;
        virtual void draw() const = 0;
    };

    ////////////////////////////////////////////////////////////////////////////
    class GameObject : public IGameObject {
        Vec2 pos;
        Vec2 velocity;
        char name[32];
        Model* model;
        // ... other members ...
        float someAccumulator;

    public:
        GameObject(float px, float py, float vx, float vy)
            : pos{ px, py }
            , velocity{ vx, vy }
            , name{ "dot object" }
            , model{ nullptr }
            , someAccumulator{ 0.0f }
        {
        }

        void advance(RotMatrix velocityRotation)
        {
            Vec2 actualVelocity = rotate(velocityRotation, velocity);
            pos.x += actualVelocity.x;
            pos.y += actualVelocity.y;
        }

        void draw() const
        {
            if (pos.x < g_viewSize && pos.y < g_viewSize)
                drawPoint(pos);
        }
    };

    ////////////////////////////////////////////////////////////////////////////
    class World {
        vector<IGameObject*> objects;

    public:
        void generateRandom()
        {
            objects.reserve(g_numPoints);
            for (int i = 0; i < g_numPoints; i++) {
                // Create an object with random position and velocity
                double px = randBetween(0.0f, g_worldSize);
                double py = randBetween(0.0f, g_worldSize);
                double vx = randBetween(-g_maxVelocity, g_maxVelocity);
                double vy = randBetween(-g_maxVelocity, g_maxVelocity);
                IGameObject* obj = nullptr;
                if ( time(NULL) > 100 )
                    obj = new GameObject(px, py, vx, vy);
                objects.push_back(obj);
            }
        }
        void advance(RotMatrix velocityRotation)
        {
            for (int i = 0; i < objects.size(); i++)
                objects[i]->advance(velocityRotation);
        }
        void draw() const
        {
            for (int i = 0; i < g_numPoints; i++)
                objects[i]->draw();
        }
    };
}

namespace ClassicOOPFlat
{
    ////////////////////////////////////////////////////////////////////////////
    class GameObject {
        Vec2 pos;
        Vec2 velocity;
        char name[32];
        Model* model;
        // ... other members ...
        float someAccumulator;

    public:
        GameObject(float px, float py, float vx, float vy)
            : pos{ px, py }
            , velocity{ vx, vy }
            , name{ "unknown obj" }
            , model{ nullptr }
            , someAccumulator{ 0.0f }
        {
        }

        void advance(RotMatrix velocityRotation)
        {
            Vec2 actualVelocity = rotate(velocityRotation, velocity);
            pos.x += actualVelocity.x;
            pos.y += actualVelocity.y;
        }

        void draw() const
        {
            if (pos.x < g_viewSize && pos.y < g_viewSize)
                drawPoint(pos);
        }
    };

    ////////////////////////////////////////////////////////////////////////////
    class World {
        vector<GameObject> objects;

    public:
        void generateRandom()
        {
            objects.reserve(g_numPoints);
            for (int i = 0; i < g_numPoints; i++) {
                // Create an object with random position and velocity
                double px = randBetween(0.0f, g_worldSize);
                double py = randBetween(0.0f, g_worldSize);
                double vx = randBetween(-g_maxVelocity, g_maxVelocity);
                double vy = randBetween(-g_maxVelocity, g_maxVelocity);
                objects.emplace_back(GameObject(px, py, vx, vy));
            }
        }
        void advance(RotMatrix velocityRotation)
        {
            for (int i = 0; i < objects.size(); i++)
                objects[i].advance(velocityRotation);
        }
        void draw() const
        {
            for (int i = 0; i < g_numPoints; i++)
                objects[i].draw();
        }
    };
}

namespace SplitStruct
{
    ////////////////////////////////////////////////////////////////////////////
    class OtherGameObjectData {
        char name[32];
        Model* model;
        // ... other members ...
        float someAccumulator;

    public:
        OtherGameObjectData()
            : name{ "unknown obj" }
            , model{ nullptr }
            , someAccumulator{ 0.0f }
        {
        }
    };

    ////////////////////////////////////////////////////////////////////////////
    class World {
        vector<Vec2> positions;
        vector<Vec2> velocities;
        vector<OtherGameObjectData> otherData;

    public:
        void generateRandom()
        {
            positions.reserve(g_numPoints);
            velocities.reserve(g_numPoints);
            otherData.reserve(g_numPoints);
            for (int i = 0; i < g_numPoints; i++) {
                // Create an object with random position and velocity
                double px = randBetween(0.0f, g_worldSize);
                double py = randBetween(0.0f, g_worldSize);
                double vx = randBetween(-g_maxVelocity, g_maxVelocity);
                double vy = randBetween(-g_maxVelocity, g_maxVelocity);
                positions.emplace_back(Vec2(px, py));
                velocities.emplace_back(Vec2(vx, vy));
                otherData.emplace_back(OtherGameObjectData());
            }
        }
        void advance(RotMatrix velocityRotation)
        {
            for (int i = 0; i < g_numPoints; i++)
            {
                Vec2 actualVelocity = rotate(velocityRotation, velocities[i]);
                positions[i].x += actualVelocity.x;
                positions[i].y += actualVelocity.y;
            }
        }
        void draw() const
        {
            for (int i = 0; i < g_numPoints; i++)
                drawPoint(positions[i]);
        }
    };
}

namespace SplitStructAndObj
{
    ////////////////////////////////////////////////////////////////////////////
    class OtherGameObjectData {
        char name[32];
        Model* model;
        // ... other members ...
        float someAccumulator;

    public:
        OtherGameObjectData()
            : name{ "unknown obj" }
            , model{ nullptr }
            , someAccumulator{ 0.0f }
        {
        }
    };

    ////////////////////////////////////////////////////////////////////////////
    class World {
        vector<Vec2> positions;
        vector<Vec2> velocities;
        vector<OtherGameObjectData> otherData;
        int endOfNear = 0;
        int updateSkipCount = 0;
        RotMatrix rotationMatSum{0,0};

    public:
        void generateRandom()
        {
            positions.reserve(g_numPoints);
            velocities.reserve(g_numPoints);
            otherData.reserve(g_numPoints);
            for (int i = 0; i < g_numPoints; i++) {
                // Create an object with random position and velocity
                double px = randBetween(0.0f, g_worldSize);
                double py = randBetween(0.0f, g_worldSize);
                double vx = randBetween(-g_maxVelocity, g_maxVelocity);
                double vy = randBetween(-g_maxVelocity, g_maxVelocity);
                // if ( px < 2*g_viewSize && py < 2*g_viewSize )
                positions.emplace_back(Vec2(px, py));
                velocities.emplace_back(Vec2(vx, vy));
                otherData.emplace_back(OtherGameObjectData());
            }
            // Bring data close to the view in front
            arrangeData();
        }
        void advance(RotMatrix velocityRotation)
        {
            // Advance only the near-view objects
            for (int i = 0; i < endOfNear; i++)
            {
                Vec2 actualVelocity = rotate(velocityRotation, velocities[i]);
                positions[i].x += actualVelocity.x;
                positions[i].y += actualVelocity.y;
            }
            // Rarely, also update the objects far away from the screen
            // pos_final = pos_0 + v*rot_1 + v*rot_2 + v*rot_3 + ...
            // pos_final = pos_0 + v*(rot1_ + rot_2 + ...)
            // pos_final = pos_0 + v*rotationMatSum
            static const int farUpdateFrequency = 100;
            rotationMatSum.cos_theta += velocityRotation.cos_theta;
            rotationMatSum.sin_theta += velocityRotation.sin_theta;
            if ( ++updateSkipCount == farUpdateFrequency )
            {
                // Walk over the data and update it; same traditional formula is used
                for (int i = endOfNear; i < g_numPoints; i++)
                {
                    Vec2 actualVelocity = rotate(rotationMatSum, velocities[i]);
                    positions[i].x += actualVelocity.x;
                    positions[i].y += actualVelocity.y;
                }
                // Start over
                rotationMatSum.cos_theta = 0.0f;
                rotationMatSum.sin_theta = 0.0f;
            }
        }
        void draw() const
        {
            for (int i = 0; i < endOfNear; i++)
                drawPoint(positions[i]);
        }

    private:
        //! Do a partial sort: bring all the elements close the view in front; push all the others
        //! to the end. Make endOfNear indicate the boundary between them
        void arrangeData()
        {
            endOfNear = 0;
            for ( int i=0; i<g_numPoints; i++ )
            {
                if (isNearView(positions[i]))
                {
                    int destIdx = endOfNear++;
                    swap(positions[i], positions[destIdx]);
                    swap(velocities[i], velocities[destIdx]);
                    swap(otherData[i], otherData[destIdx]);
                }
            }
        }

        bool isNearView(Vec2 pos)
        {
            return pos.x < 4*g_viewSize && pos.y < 4*g_viewSize;
        }
    };
}

typedef ClassicOOP::World WorldType;        // 180+ ms
// typedef ClassicOOPFlat::World WorldType;    // 140+ ms
// typedef SplitStruct::World WorldType;       // 30+ ms
// typedef SplitStructAndObj::World WorldType; // < 16 ms

WorldType g_world;

void displayOneFrame()
{
    // First advance the world
    // Always rotate the velocities, for a nice effect
    static float rotAngle = 0.0f;
    rotAngle += 0.01;
    g_world.advance(getRotMatrix(rotAngle)); //  try commenting it out

    // Now draw the points

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    glPointSize(4.0f);

    glBegin(GL_POINTS);
    glColor3f(1.0f, 0.9f, 0.9f);

    // The actual draw - try commenting it out
    g_world.draw();

    glEnd();
    glFlush();

    // Calculate the frames per second
    static int frame = 0;
    static int curTime = 0;
    static int lastSecTime = 0;
    frame++;
    curTime = glutGet(GLUT_ELAPSED_TIME);
    if (curTime - lastSecTime > 1000) {
        static char title[256];
        float fps = frame * 1000.0 / (curTime - lastSecTime);
        int frameTime = int(1000.0f / fps);
        sprintf(title, "GameObject Test - FPS: %4.2f - avg frameTime: %d ms", fps, frameTime);
        glutSetWindowTitle(title);
        lastSecTime = curTime;
        frame = 0;
    }
}

int main(int argc, char** argv)
{
    srand(1);

    printf("sizeof(GameObject)=%d\n", (int) sizeof(ClassicOOP::GameObject));
    printf("sizeof(Vec2)=%d\n", (int) sizeof(Vec2));

    // Initialize the world with random data
    g_world.generateRandom();

    glutInit(&argc, argv);
    glutInitWindowSize(g_viewSize, g_viewSize);
    glutCreateWindow("GameObject Test");
    glutDisplayFunc(displayOneFrame);
    glutIdleFunc(displayOneFrame);
    glutMainLoop();
    return 0;
}
