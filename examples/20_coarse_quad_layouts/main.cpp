/* This sample program computes the coarse quad
 * decomposition of a given quadrilateral mesh.
 *
 * Enjoy!
*/

#include <QApplication>
#include <QSlider>
#include <QPushButton>
#include <QGridLayout>
#include <cinolib/meshes/meshes.h>
#include <cinolib/profiler.h>
#include <cinolib/gui/qt/qt_gui_tools.h>
#include <cinolib/coarse_layout.h>
#include <cinolib/drawable_sphere.h>

//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
using namespace cinolib;
using namespace std;

int main(int argc, char **argv)
{
    QApplication a(argc, argv);

    string s = (argc==2) ? string(argv[1]) : string(DATA_PATH) + "cubespikes.obj";
    DrawableQuadmesh<> m(s.c_str());

    Profiler profiler;
    profiler.push("coarse layout");
    compute_coarse_quad_layout(m);
    profiler.pop();
    m.poly_color_wrt_label();
    m.show_marked_edge_color(Color::BLACK());
    m.show_marked_edge_width(3);
    m.edge_set_alpha(0.5);

    GLcanvas gui;
    gui.push_obj(&m);
    gui.show();

    for(uint vid=0; vid<m.num_verts(); ++vid)
        if(m.vert_is_singular(vid)) gui.push_obj(new DrawableSphere(m.vert(vid), 0.5));
    
    // CMD+1 to show mesh controls.
    SurfaceMeshControlPanel<DrawableQuadmesh<>> panel(&m, &gui);
    QApplication::connect(new QShortcut(QKeySequence(Qt::CTRL+Qt::Key_1), &gui), &QShortcut::activated, [&](){panel.show();});

    return a.exec();
}
