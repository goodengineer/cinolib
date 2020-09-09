/* This is a base application for cinolib (https://github.com/maxicino/cinolib).
 * It will show a GL canvas (and associated control panel) to interact
 * with a tetrahedral mesh.
 * In case you don't need a GUI, you can drop the "Drawable" prefix from the mesh data type.
 * What you will get is a lighter yet fully operational mesh data structure, just
 * without the burden of OpenGL code necessary for rendering!
 *
 * Enjoy!
*/

#include <QApplication>
#include <cinolib/meshes/meshes.h>
#include <cinolib/gui/qt/qt_gui_tools.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

using namespace cinolib;
using namespace std;

int main(int argc, char **argv)
{
    QApplication a(argc, argv);

    string s = (argc==2) ? string(argv[1]) : string(DATA_PATH) + "/sphere.mesh";

    DrawableTetmesh<> m(s.c_str());

    GLcanvas gui;
    gui.push_obj(&m);
    gui.show();

    // CMD+1 to show mesh controls.
    VolumeMeshControlPanel<DrawableTetmesh<>> panel(&m, &gui);
    QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_1), &gui), &QShortcut::activated, [&](){panel.show();});

    return a.exec();
}
