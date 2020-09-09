/* This is a base application for cinolib (https://github.com/maxicino/cinolib).
 * It will show a GL canvas (and associated control panel) to interact
 * with a quadrilateral mesh.
 *
 * Enjoy!
*/

#include <QApplication>
#include <cinolib/meshes/meshes.h>
#include <cinolib/gui/qt/qt_gui_tools.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#define DRAW_MESH TRUE

int main(int argc, char **argv)
{
    using namespace cinolib;

    QApplication a(argc, argv);

    std::string s = (argc==2) ? std::string(argv[1]) : std::string(DATA_PATH) + "/cubespikes.obj";

    #ifdef DRAW_MESH
    DrawableQuadmesh<> m(s.c_str());

    GLcanvas gui;
    gui.push_obj(&m);
    gui.show();

    // CMD+1 to show mesh controls.
    SurfaceMeshControlPanel<DrawableQuadmesh<>> panel(&m, &gui);
    QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_1), &gui), &QShortcut::activated, [&](){panel.show();});

    #else
    Quadmesh<> m(s.c_str());
    // In case you don't need a GUI, you can drop the "Drawable" prefix from the mesh data type.
    // What you will get is a lighter yet fully operational mesh data structure, just
    // without the burden of OpenGL code necessary for rendering!
    //Your model(mesh) specific processing code goes here:
    
    #endif
    
    return a.exec();
}
