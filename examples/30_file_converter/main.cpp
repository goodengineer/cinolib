/* This is a simple file converter tool for all the file formats supported by cinolib
 *
 * Enjoy!
*/

#include <cinolib/io/read_write.h>
#include <cinolib/string_utilities.h>
#include <cinolib/stl_container_utilities.h>
#include <cinolib/meshes/trimesh.h>

using namespace cinolib;
using namespace std;
//::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

int main(int argc, char **argv)
{
    if(argc!=3)
    {
        cout<<"\n\nusage:\n\tfile_converter input output\n\n"<<endl;;
        return -1;
    }

    vector<vec3d>        verts;
    vector<vector<uint>> faces;
    vector<vector<uint>> polys;
    vector<vector<bool>> winding;

    // input extension
    string ext = get_file_extension(argv[1]);

    if(ext.compare("OBJ")==0 || ext.compare("obj")==0)
        read_OBJ(argv[1], verts, polys);
    else if(ext.compare("OFF")==0 || ext.compare("off")==0)
        read_OFF(argv[1], verts, polys);
    else if(ext.compare("STL")==0 || ext.compare("stl")==0)
    {
        vector<uint> tris;
        read_STL(argv[1], verts, tris);
        polys = polys_from_serialized_vids(tris,3);
    }
    else if(ext.compare("HEDRA")==0 || ext.compare("hedra")==0)
        read_HEDRA(argv[1], verts, faces, polys, winding);
    else if(ext.compare("HYBRID")==0 || ext.compare("hybrid")==0)
        read_HYBDRID(argv[1], verts, faces, polys, winding);
    else if(ext.compare("MESH")==0 || ext.compare("mesh")==0)
        read_MESH(argv[1], verts, polys);
    else if(ext.compare("TET")==0 || ext.compare("tet")==0)
        read_TET(argv[1], verts, polys);
    else if(ext.compare("VTU")==0 || ext.compare("vtu")==0)
        read_VTU(argv[1], verts, polys);
    else if(ext.compare("VTK")==0 || ext.compare("vtk")==0)
        read_VTK(argv[1], verts, polys);
    else
    {
        cout<<"ERROR: unknown input format"<<endl;;
        return -1;
    }

    // output extension
    ext = get_file_extension(argv[2]);

    if(ext.compare("OBJ")==0 || ext.compare("obj")==0)
        write_OBJ(argv[2], serialized_xyz_from_vec3d(verts), polys);
    else if(ext.compare("OFF")==0 || ext.compare("off")==0)
        write_OFF(argv[2], serialized_xyz_from_vec3d(verts), polys);
    else if(ext.compare("STL")==0 || ext.compare("stl")==0)
    {
        Trimesh<> tmp(verts,polys);
        tmp.save(argv[2]); // STL needs surface normals, so I make the mesh
    }
    else if(ext.compare("NODE")==0 || ext.compare("node")==0 ||
            ext.compare("ELE")==0  || ext.compare("ele")==0)
    {
        string tmp = get_file_name(argv[2], false);
        write_NODE_ELE(tmp.c_str(), verts, polys);
    }
    else if(ext.compare("HEDRA")==0 || ext.compare("hedra")==0)
        write_HEDRA(argv[2], verts, faces, polys, winding);
    else if(ext.compare("MESH")==0 || ext.compare("mesh")==0)
        write_MESH(argv[2], verts, polys);
    else if(ext.compare("TET")==0 || ext.compare("tet")==0)
        write_TET(argv[2], verts, polys);
    else if(ext.compare("VTU")==0 || ext.compare("vtu")==0)
        write_VTU(argv[2], verts, polys);
    else if(ext.compare("VTK")==0 || ext.compare("vtk")==0)
        write_VTK(argv[2], verts, polys);
    else
    {
        cout<<"ERROR: unknown output format"<<endl;;
        return -1;
    }

    return 0;
}
