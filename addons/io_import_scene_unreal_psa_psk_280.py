# ##### BEGIN GPL LICENSE BLOCK #####
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software Foundation,
#  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ##### END GPL LICENSE BLOCK #####

bl_info = {
    "name": "Import Unreal Skeleton Mesh (.psk)/Animation Set (.psa) (280)",
    "author": "Darknet, flufy3d, camg188, befzz",
    "version": (2, 8, 0),
    "blender": (2, 80, 0),
    "location": "File > Import > Skeleton Mesh (.psk)/Animation Set (.psa) OR View3D > Tool Shelf (key T) > Misc. tab",
    "description": "Import Skeleton Mesh / Animation Data",
    "warning": "",
    "wiki_url": "https://github.com/Befzz/blender3d_import_psk_psa",
    "category": "Import-Export",
    "tracker_url": "https://github.com/Befzz/blender3d_import_psk_psa/issues"
}

"""
Version': '2.0' ported by Darknet

Unreal Tournament PSK file to Blender mesh converter V1.0
Author: D.M. Sturgeon (camg188 at the elYsium forum), ported by Darknet
Imports a *psk file to a new mesh

-No UV Texutre
-No Weight
-No Armature Bones
-No Material ID
-Export Text Log From Current Location File (Bool )
"""

"""
Version': '2.7.*' edited by befzz
Github: https://github.com/Befzz/blender3d_import_psk_psa
- Pskx support
- Animation import updated (bone orientation now works)
- Skeleton import: auto-size, auto-orient bones
- UVmap, mesh, weights, etc. import revised
- Extra UVs import

- No Scale support. (no test material)
- No smoothing groups (not exported by umodel)
"""

"""
Version': '2.8.0' edited by floxay
- Vertex normals import (VTXNORMS chunk)
        (requires custom UEViewer build /at the moment/)
"""

# https://github.com/gildor2/UModel/blob/master/Exporters/Psk.h

import bpy
import re
from mathutils import Vector, Matrix, Quaternion
from bpy.props import (FloatProperty,
                        StringProperty,
                        BoolProperty,
                        EnumProperty,
                        PointerProperty )
                      
from struct import unpack, unpack_from, Struct
import time

#DEV
# from mathutils import *
# from math import *


def util_obj_link(context, obj):
    # return bpy.context.scene_collection.objects.link(obj)
    # bpy.context.view_layer.collections[0].collection.objects.link(obj)
    # return bpy.context.collection.objects.link(obj)
    # bpy.data.scenes[0].collection.objects.link(obj)
    context.collection.objects.link(obj)

def util_obj_select(context, obj, action = 'SELECT'):
    # if obj.name in bpy.data.scenes[0].view_layers[0].objects:
    if obj.name in context.view_layer.objects:
        return obj.select_set(action == 'SELECT')
    else:
        print('Warning: util_obj_select: Object not in "context.view_layer.objects"')

def util_obj_set_active(context, obj):
    # bpy.context.view_layer.objects.active = obj
    # bpy.data.scenes[0].view_layers[0].objects.active = obj
    context.view_layer.objects.active = obj

def util_get_scene(context):
    return context.scene

def get_uv_layers(mesh_obj):
    return mesh_obj.uv_layers
    
def obj_select_get(obj):
    return obj.select_get()


def utils_set_mode(mode):
    if bpy.ops.object.mode_set.poll():
        bpy.ops.object.mode_set(mode = mode, toggle = False)
    # else:
        # bpy.ops.object.mode_set(mode = mode, toggle = False)
        #dev

# since names have type ANSICHAR(signed char) - using cp1251(or 'ASCII'?)
def util_bytes_to_str(in_bytes):
    return in_bytes.rstrip(b'\x00').decode(encoding = 'cp1252', errors = 'replace')

class class_psk_bone:
    name = ""
    
    parent = None
    
    bone_index = 0
    parent_index = 0
    
    # scale = []
    
    mat_world = None
    mat_world_rot = None
    
    orig_quat = None
    orig_loc = None
    
    children = None
    
    have_weight_data = False

# TODO simplify?
def util_select_all(select):
    if select:
        actionString = 'SELECT'
    else:
        actionString = 'DESELECT'

    if bpy.ops.object.select_all.poll():
        bpy.ops.object.select_all(action = actionString)

    if bpy.ops.mesh.select_all.poll():
        bpy.ops.mesh.select_all(action = actionString)

    if bpy.ops.pose.select_all.poll():
        bpy.ops.pose.select_all(action = actionString)

        
def util_ui_show_msg(msg):
    bpy.ops.pskpsa.message('INVOKE_DEFAULT', message = msg)
        
        
PSKPSA_FILE_HEADER = {
    'psk':b'ACTRHEAD\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00',
    'psa':b'ANIMHEAD\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'
}


def util_is_header_valid(filename, file_ext, chunk_id, error_callback):
    '''Return True if chunk_id is a valid psk/psa (file_ext) 'magick number'.'''
    if chunk_id != PSKPSA_FILE_HEADER[file_ext]:
        error_callback(
            "File %s is not a %s file. (header mismach)\nExpected: %s \nPresent %s"  % ( 
                filename, file_ext,
                PSKPSA_FILE_HEADER[file_ext], chunk_id)
        )    
        return False
    return True
    
    
def util_gen_name_part(filepath):
    '''Return file name without extension'''
    return re.match(r'.*[/\\]([^/\\]+?)(\..{2,5})?$', filepath).group(1)
                
                
def vec_to_axis_vec(vec_in, vec_out):
    '''Make **vec_out** to be an axis-aligned unit vector that is closest to vec_in. (basis?)'''
    x, y, z = vec_in
    if abs(x) > abs(y):
        if abs(x) > abs(z):
            vec_out.x = 1 if x >= 0 else -1
        else:
            vec_out.z = 1 if z >= 0 else -1
    else:
        if abs(y) > abs(z):
            vec_out.y = 1 if y >= 0 else -1
        else:
            vec_out.z = 1 if z >= 0 else -1
            
            
def calc_bone_rotation(psk_bone, bone_len, bDirectly, avg_bone_len):
    children = psk_bone.children

    vecy = Vector((0.0, 1.0, 0.0))
    quat = Quaternion((1.0, 0.0, 0.0, 0.0))
    axis_vec = Vector()
    
    # bone with 0 children (orphan bone)
    if len(children) == 0:
        # Single bone. ALONE.
        if psk_bone.parent == None:
            return (bone_len, quat)
            
        elif bDirectly:
                        # @
            # axis_vec = psk_bone.orig_quat * psk_bone.orig_loc
            axis_vec = psk_bone.orig_loc.copy()
            axis_vec.rotate( psk_bone.orig_quat )

        else:
            # bone Head near parent Head?
            if psk_bone.orig_loc.length < 0.1 * avg_bone_len:
                # @
                # vec_to_axis_vec(psk_bone.orig_quat.conjugated() * psk_bone.parent.axis_vec, axis_vec)
                v = psk_bone.parent.axis_vec.copy()
                v.rotate( psk_bone.orig_quat.conjugated() )
                vec_to_axis_vec(v, axis_vec)

                # reorient bone to other axis bychanging our base Y vec...
                # this is not tested well
                vecy = Vector((1.0, 0.0, 0.0))
            else:
                # @
                # vec_to_axis_vec(psk_bone.orig_quat.conjugated() * psk_bone.parent.axis_vec, axis_vec)
                v = psk_bone.parent.axis_vec.copy()
                v.rotate( psk_bone.orig_quat.conjugated() )
                vec_to_axis_vec(v, axis_vec)

        return (bone_len, vecy.rotation_difference(axis_vec))
        
    # bone with > 0 children BUT only 1 non orphan bone ( reorient to it! )
    if bDirectly and len(children) > 1:
    
        childs_with_childs = 0
                    
        for child in filter(lambda c: len(c.children), children):
 
            childs_with_childs += 1
            
            if childs_with_childs > 1:
                break
                
            candidate = child
                    
        if childs_with_childs == 1:
            # print('candidate',psk_bone.name,candidate.name)
            return (len(candidate.orig_loc), vecy.rotation_difference(candidate.orig_loc))
            
    # bone with > 0 children
    sumvec = Vector()
    sumlen = 0
    
    for child in children:
        sumvec += (child.orig_loc)
        sumlen += child.orig_loc.length
    sumlen /= len(children)
    sumlen = max(sumlen, 0.01)
    
    if bDirectly:
        return (sumlen, vecy.rotation_difference(sumvec))
    
    vec_to_axis_vec(sumvec, axis_vec)
    psk_bone.axis_vec = axis_vec
    return (sumlen, vecy.rotation_difference(axis_vec))

    
def __pass(*args,**kwargs):
    pass


def util_check_file_header(file, ftype):
    header_bytes = file.read(32)
    
    if len(header_bytes) < 32:
        return False
        
    if not header_bytes.startswith( PSKPSA_FILE_HEADER[ftype] ):
        return False
        
    return True
        
    
def color_linear_to_srgb(c):
    """
    Convert from linear to sRGB color space.
    Source: Cycles addon implementation, node_color.h.
    """
    if c < 0.0031308:
        return 0.0 if c < 0.0 else c * 12.92
    else:
        return 1.055 * pow(c, 1.0 / 2.4) - 0.055
        
def pskimport(filepath,
        context = None,
        bImportmesh = True,
        bImportbone = True,
        bSpltiUVdata = False,
        fBonesize = 5.0,
        fBonesizeRatio = 0.6,
        bDontInvertRoot = True,
        bReorientBones = False,
        bReorientDirectly = False,
        bScaleDown = True,
        bToSRGB = True,
        error_callback = None):
    '''
    Import mesh and skeleton from .psk/.pskx files
    
    Args:
        bReorientBones:
            Axis based bone orientation to children
            
        error_callback:
            Called when importing is failed.
            
            error_callback = lambda msg: print('reason:', msg)
            
    '''
    if not hasattr( error_callback, '__call__'):
        # error_callback = __pass
        error_callback = print
        
    # ref_time = time.process_time()
    if not bImportbone and not bImportmesh:
        error_callback("Nothing to do.\nSet something for import.")
        return False
        
    print ("-----------------------------------------------")
    print ("---------EXECUTING PSK PYTHON IMPORTER---------")
    print ("-----------------------------------------------")

    #file may not exist
    try:
        file = open(filepath,'rb')
    except IOError:
        error_callback('Error while opening file for reading:\n  "'+filepath+'"')
        return False

    if not util_check_file_header(file, 'psk'):
        error_callback('Not psk file:\n  "'+filepath+'"')
        return False
        
    Vertices = None
    Wedges = None
    Faces = None
    UV_by_face = None
    Materials = None
    Bones = None
    Weights = None
    VertexColors = None
    Extrauvs = []
    Normals = None
    WedgeIdx_by_faceIdx = None
     
    if not context:
        context = bpy.context
    #================================================================================================== 
    # Materials   MaterialNameRaw | TextureIndex | PolyFlags | AuxMaterial | AuxFlags |  LodBias | LodStyle 
    # Only Name is usable.
    def read_materials():
    
        nonlocal Materials
        
        Materials = []
        
        for counter in range(chunk_datacount):

            (MaterialNameRaw,) = unpack_from('64s24x', chunk_data, chunk_datasize * counter)
            
            Materials.append( util_bytes_to_str( MaterialNameRaw ) )
            
            
    #================================================================================================== 
    # Faces WdgIdx1 | WdgIdx2 | WdgIdx3 | MatIdx | AuxMatIdx | SmthGrp
    def read_faces():
        
        if not bImportmesh:
            return True
        
        nonlocal Faces, UV_by_face, WedgeIdx_by_faceIdx

        UV_by_face = [None] * chunk_datacount
        Faces = [None] * chunk_datacount
        WedgeIdx_by_faceIdx = [None] * chunk_datacount
        
        if len(Wedges) > 65536:
            unpack_format = '=IIIBBI'
        else:
            unpack_format = '=HHHBBI'
            
        unpack_data = Struct(unpack_format).unpack_from
        
        for counter in range(chunk_datacount):
            (WdgIdx1, WdgIdx2, WdgIdx3,
             MatIndex, 
             AuxMatIndex, #unused
             SmoothingGroup # Umodel is not exporting SmoothingGroups
             ) = unpack_data(chunk_data, counter * chunk_datasize)
             
            # looks ugly
            # Wedges is (point_index, u, v, MatIdx)
            ((vertid0, u0, v0, matid0), (vertid1, u1, v1, matid1), (vertid2, u2, v2, matid2)) = Wedges[WdgIdx1], Wedges[WdgIdx2], Wedges[WdgIdx3]
            
            # note order: C,B,A
            # Faces[counter] = (vertid2,  vertid1, vertid0)

            Faces[counter] = (vertid1,  vertid0, vertid2)
            # Faces[counter] = (vertid1,  vertid2, vertid0)
            # Faces[counter] = (vertid0,  vertid1, vertid2)
            
            # uv = ( ( u2, 1.0 - v2 ), ( u1, 1.0 - v1 ), ( u0, 1.0 - v0 ) )
            uv = ( ( u1, 1.0 - v1 ), ( u0, 1.0 - v0 ), ( u2, 1.0 - v2 ) )
            
            # Mapping: FaceIndex <=> UV data <=> FaceMatIndex
            UV_by_face[counter] = (uv, MatIndex, (matid2, matid1, matid0))
            
            # We need this for EXTRA UVs
            WedgeIdx_by_faceIdx[counter] = (WdgIdx3, WdgIdx2, WdgIdx1)

            
    #==================================================================================================
    # Vertices X | Y | Z
    def read_vertices():
        
        if not bImportmesh:
            return True
            
        nonlocal Vertices
        
        Vertices = [None] * chunk_datacount
        
        unpack_data = Struct('3f').unpack_from
        
        if bScaleDown:
            for counter in range( chunk_datacount ):
                (vec_x, vec_y, vec_z) = unpack_data(chunk_data, counter * chunk_datasize)
                Vertices[counter]  = (vec_x*0.01, vec_y*0.01, vec_z*0.01)
                # equal to gltf
                # Vertices[counter]  = (vec_x*0.01, vec_z*0.01, -vec_y*0.01)
        else:
            for counter in range( chunk_datacount ):
                Vertices[counter]  =  unpack_data(chunk_data, counter * chunk_datasize)
            
            
    #================================================================================================== 
    # Wedges (UV)   VertexId |  U |  V | MatIdx 
    def read_wedges():
    
        if not bImportmesh:
            return True
            
        nonlocal Wedges
        
        Wedges = [None] * chunk_datacount
        
        unpack_data = Struct('=IffBxxx').unpack_from
        
        for counter in range( chunk_datacount ):
            (vertex_id,
             u, v,
             material_index) = unpack_data( chunk_data, counter * chunk_datasize )
             
            # print(vertex_id, u, v, material_index)
            # Wedges[counter] = (vertex_id, u, v, material_index)
            Wedges[counter] = [vertex_id, u, v, material_index]
            
    #================================================================================================== 
    # Bones (VBone .. VJointPos ) Name|Flgs|NumChld|PrntIdx|Qw|Qx|Qy|Qz|LocX|LocY|LocZ|Lngth|XSize|YSize|ZSize
    def read_bones():
    
        nonlocal Bones, bImportbone
        
        if chunk_datacount == 0:
            bImportbone = False
            
        if bImportbone:
            # unpack_data = Struct('64s3i11f').unpack_from
            unpack_data = Struct('64s3i7f16x').unpack_from
        else:
            unpack_data = Struct('64s56x').unpack_from
            
        Bones = [None] * chunk_datacount
        
        for counter in range( chunk_datacount ):
            Bones[counter] = unpack_data( chunk_data, chunk_datasize * counter)
          
    
    #================================================================================================== 
    # Influences (Bone Weight) (VRawBoneInfluence) ( Weight | PntIdx | BoneIdx)
    def read_weights():

        nonlocal Weights
        
        if not bImportmesh:
            return True
            
        Weights = [None] * chunk_datacount
        
        unpack_data = Struct('fii').unpack_from
        
        for counter in range(chunk_datacount):
            Weights[counter] = unpack_data(chunk_data, chunk_datasize * counter)
             
    #================================================================================================== 
    # Vertex colors. R G B A bytes. NOTE: it is Wedge color.(uses Wedges index)
    def read_vertex_colors():
    
        nonlocal VertexColors
        
        unpack_data = Struct("=4B").unpack_from
        
        VertexColors = [None] * chunk_datacount
        
        for counter in range( chunk_datacount ):
            VertexColors[counter] = unpack_data(chunk_data, chunk_datasize * counter) 
            
    
    #================================================================================================== 
    # Extra UV. U | V
    def read_extrauvs():

        unpack_data = Struct("=2f").unpack_from
        
        uvdata = [None] * chunk_datacount
        
        for counter in range( chunk_datacount ):
            uvdata[counter] = unpack_data(chunk_data, chunk_datasize * counter) 
            
        Extrauvs.append(uvdata)

    #==================================================================================================
    # Vertex Normals NX | NY | NZ
    def read_normals():
        if not bImportmesh:
            return True

        nonlocal Normals
        Normals = [None] * chunk_datacount

        unpack_data = Struct('3f').unpack_from

        for counter in range(chunk_datacount):
            Normals[counter] = unpack_data(chunk_data, counter * chunk_datasize)
 
             
    CHUNKS_HANDLERS = {
        'PNTS0000': read_vertices,
        'VTXW0000': read_wedges,
        'VTXW3200': read_wedges,#?
        'FACE0000': read_faces,
        'FACE3200': read_faces,
        'MATT0000': read_materials,
        'REFSKELT': read_bones,
        'REFSKEL0': read_bones, #?
        'RAWW0000': read_weights,
        'RAWWEIGH': read_weights,
        'VERTEXCO': read_vertex_colors, # VERTEXCOLOR
        'EXTRAUVS': read_extrauvs,
        'VTXNORMS': read_normals
    }
    
    #===================================================================================================
    # File. Read all needed data.
    #         VChunkHeader Struct
    # ChunkID|TypeFlag|DataSize|DataCount
    # 0      |1       |2       |3
    
    while True:
    
        header_bytes = file.read(32)
        
        if len(header_bytes) < 32:
            
            if len(header_bytes) != 0:
                error_callback("Unexpected end of file.(%s/32 bytes)" % len(header_bytes))
            break
            
        (chunk_id, chunk_type, chunk_datasize, chunk_datacount) = unpack('20s3i', header_bytes)
        
        chunk_id_str = util_bytes_to_str(chunk_id)
        chunk_id_str = chunk_id_str[:8]
        
        if chunk_id_str in CHUNKS_HANDLERS:
        
            chunk_data = file.read( chunk_datasize * chunk_datacount)
            
            if len(chunk_data) < chunk_datasize * chunk_datacount:
                error_callback('Psk chunk %s is broken.' % chunk_id_str)
                return False
                
            CHUNKS_HANDLERS[chunk_id_str]()
            
        else:
        
            print('Unknown chunk: ', chunk_id_str)
            file.seek(chunk_datasize * chunk_datacount, 1)
            
            
        # print(chunk_id_str, chunk_datacount)
            
    file.close()
         
    print(" Importing file:", filepath)
    
    if not bImportmesh and (Bones is None or len(Bones) == 0):
        error_callback("Psk: no skeleton data.")
        return False

    MAX_UVS = 8
    NAME_UV_PREFIX = "UV"
    
    # file name w/out extension
    gen_name_part = util_gen_name_part(filepath)
    gen_names = {
        'armature_object':  gen_name_part + '.ao',
        'armature_data':    gen_name_part + '.ad',
            'mesh_object':  gen_name_part + '.mo',
            'mesh_data':    gen_name_part + '.md'
    }
    
    if bImportmesh:
        mesh_data = bpy.data.meshes.new(gen_names['mesh_data'])
        mesh_obj = bpy.data.objects.new(gen_names['mesh_object'], mesh_data)
        
    
    #==================================================================================================
    # UV. Prepare
    if bImportmesh:
        if bSpltiUVdata:
        # store how much each "matrial index" have vertices
        
            uv_mat_ids = {}
            
            for (_, _, _, material_index) in Wedges:
            
                if not (material_index in uv_mat_ids):
                    uv_mat_ids[material_index] = 1
                else:
                    uv_mat_ids[material_index] += 1
                    
                    
            # if we have more UV material indexes than blender UV maps, then...
            if bSpltiUVdata and len(uv_mat_ids) > MAX_UVS :
            
                uv_mat_ids_len = len(uv_mat_ids)
                
                print('UVs: %s out of %s is combined in a first UV map(%s0)' % (uv_mat_ids_len - 8, uv_mat_ids_len, NAME_UV_PREFIX))
                
                mat_idx_proxy = [0] * len(uv_mat_ids)
                
                counts_sorted = sorted(uv_mat_ids.values(), reverse = True)
                
                new_mat_index = MAX_UVS - 1
                
                for c in counts_sorted:
                    for mat_idx, counts in uv_mat_ids.items():
                        if c == counts:
                            mat_idx_proxy[mat_idx] = new_mat_index
                            if new_mat_index > 0:
                                new_mat_index -= 1
                            # print('MatIdx remap: %s > %s' % (mat_idx,new_mat_index))
                            
                for i in range(len(Wedges)):
                    Wedges[i][3] = mat_idx_proxy[Wedges[i][3]]

        # print('Wedges:', chunk_datacount)
        # print('uv_mat_ids', uv_mat_ids)
        # print('uv_mat_ids', uv_mat_ids)
        # for w in Wedges:
        
    if bImportmesh:
        # print("-- Materials -- (index, name, faces)")
        blen_materials = []
        for materialname in Materials:
            matdata = bpy.data.materials.get(materialname)
            
            if matdata is None:
                matdata = bpy.data.materials.new( materialname )
            # matdata = bpy.data.materials.new( materialname )
                
            blen_materials.append( matdata )
            mesh_data.materials.append( matdata )
            # print(counter,materialname,TextureIndex)
            # if mat_groups.get(counter) is not None:
                # print("%i: %s" % (counter, materialname), len(mat_groups[counter]))

    #==================================================================================================
    # Prepare bone data
    def init_psk_bone(i, psk_bones, name_raw):
        psk_bone = class_psk_bone()
        psk_bone.children = []
        psk_bone.name = util_bytes_to_str(name_raw)
        psk_bones[i] = psk_bone
        return psk_bone

    psk_bone_name_toolong = False
        
    # indexed by bone index. array of psk_bone
    psk_bones = [None] * len(Bones)
    
    if not bImportbone: #data needed for mesh-only import
    
        for counter,(name_raw,) in enumerate(Bones):
            init_psk_bone(counter, psk_bones, name_raw)
            
    if bImportbone:  #else?
    
        # average bone length
        sum_bone_pos = 0
        
        for counter, (name_raw, flags, NumChildren, ParentIndex, #0 1 2 3
             quat_x, quat_y, quat_z, quat_w,            #4 5 6 7
             vec_x, vec_y, vec_z
            #  ,                       #8 9 10
            #  joint_length,                              #11
            #  scale_x, scale_y, scale_z
             ) in enumerate(Bones):
        
            psk_bone = init_psk_bone(counter, psk_bones, name_raw)
            
            psk_bone.bone_index = counter
            psk_bone.parent_index = ParentIndex
            
            # Tested. 64 is getting cut to 63
            if len(psk_bone.name) > 63:
                psk_bone_name_toolong = True
                # print('Warning. Bone name is too long:', psk_bone.name)

            # make sure we have valid parent_index
            if psk_bone.parent_index < 0:
                psk_bone.parent_index = 0
            
            # psk_bone.scale = (scale_x, scale_y, scale_z)
            # print("%s: %03f %03f | %f" % (psk_bone.name, scale_x, scale_y, joint_length),scale_x)
            # print("%s:" % (psk_bone.name), vec_x, quat_x)

            # store bind pose to make it available for psa-import via CustomProperty of the Blender bone
            psk_bone.orig_quat = Quaternion((quat_w, quat_x, quat_y, quat_z))

            if bScaleDown:
                psk_bone.orig_loc = Vector((vec_x * 0.01, vec_y * 0.01, vec_z * 0.01))
            else:
                psk_bone.orig_loc = Vector((vec_x, vec_y, vec_z))

            # root bone must have parent_index = 0 and selfindex = 0
            if psk_bone.parent_index == 0 and psk_bone.bone_index == psk_bone.parent_index:
                if bDontInvertRoot:
                    psk_bone.mat_world_rot = psk_bone.orig_quat.to_matrix()
                else:
                    psk_bone.mat_world_rot = psk_bone.orig_quat.conjugated().to_matrix()
                psk_bone.mat_world = Matrix.Translation(psk_bone.orig_loc)

            sum_bone_pos += psk_bone.orig_loc.length
            
    
    #==================================================================================================
    # Bones. Calc World-space matrix
    
        # TODO optimize math.
        for psk_bone in psk_bones:
                
            if psk_bone.parent_index == 0:
                if psk_bone.bone_index == 0:
                    psk_bone.parent = None
                    continue
                    
            parent = psk_bones[psk_bone.parent_index]
            
            psk_bone.parent = parent
            
            parent.children.append(psk_bone)
            
            # mat_world -     world space bone matrix WITHOUT own rotation
            # mat_world_rot - world space bone rotation WITH own rotation

            # psk_bone.mat_world = parent.mat_world_rot.to_4x4()
            # psk_bone.mat_world.translation = parent.mat_world.translation + parent.mat_world_rot * psk_bone.orig_loc
            # psk_bone.mat_world_rot = parent.mat_world_rot * psk_bone.orig_quat.conjugated().to_matrix()

            psk_bone.mat_world = parent.mat_world_rot.to_4x4()

            v = psk_bone.orig_loc.copy()
            v.rotate( parent.mat_world_rot )
            psk_bone.mat_world.translation = parent.mat_world.translation + v


            psk_bone.mat_world_rot = psk_bone.orig_quat.conjugated().to_matrix()
            psk_bone.mat_world_rot.rotate( parent.mat_world_rot )


            # psk_bone.mat_world =  ( parent.mat_world_rot.to_4x4() * psk_bone.trans)
            # psk_bone.mat_world.translation += parent.mat_world.translation
            # psk_bone.mat_world_rot = parent.mat_world_rot * psk_bone.orig_quat.conjugated().to_matrix()
            
    
    #==================================================================================================
    # Skeleton. Prepare.
            
        armature_data = bpy.data.armatures.new(gen_names['armature_data'])
        armature_obj = bpy.data.objects.new(gen_names['armature_object'], armature_data)
        # TODO: options for axes and x_ray?
        armature_data.show_axes = False

        armature_data.display_type = 'STICK'
        armature_obj.show_in_front = True

        util_obj_link(context, armature_obj)

        util_select_all(False)
        util_obj_select(context, armature_obj)
        util_obj_set_active(context, armature_obj)
        
        utils_set_mode('EDIT')
        
        
        sum_bone_pos /= len(Bones) # average
        sum_bone_pos *= fBonesizeRatio # corrected
        
        # bone_size_choosen = max(0.01, round((min(sum_bone_pos, fBonesize))))
        bone_size_choosen = max(0.01, round((min(sum_bone_pos, fBonesize))*100)/100)
        # bone_size_choosen = max(0.01, min(sum_bone_pos, fBonesize))
        # print("Bonesize %f | old: %f round: %f" % (bone_size_choosen, max(0.01, min(sum_bone_pos, fBonesize)),max(0.01, round((min(sum_bone_pos, fBonesize))*100)/100)))

        if not bReorientBones:
            new_bone_size = bone_size_choosen
        
    #==================================================================================================
    # Skeleton. Build.
        if psk_bone_name_toolong:
            print('Warning. Some bones will be renamed(names are too long). Animation import may be broken.')
            for psk_bone in psk_bones:

                # TODO too long name cutting options?
                orig_long_name = psk_bone.name

                # Blender will cut the name here (>63 chars)
                edit_bone = armature_obj.data.edit_bones.new(psk_bone.name)
                edit_bone["orig_long_name"] = orig_long_name

                # if orig_long_name != edit_bone.name:
                #     print('--')
                #     print(len(orig_long_name),orig_long_name)
                #     print(len(edit_bone.name),edit_bone.name)

                # Use the bone name made by blender (.001 , .002 etc.)
                psk_bone.name = edit_bone.name

        else:
            for psk_bone in psk_bones:
                edit_bone = armature_obj.data.edit_bones.new(psk_bone.name)
                psk_bone.name = edit_bone.name

        for psk_bone in psk_bones:
            edit_bone = armature_obj.data.edit_bones[psk_bone.name]

            armature_obj.data.edit_bones.active = edit_bone

            if psk_bone.parent is not None:
                edit_bone.parent = armature_obj.data.edit_bones[psk_bone.parent.name]
            else:
                if bDontInvertRoot:
                    psk_bone.orig_quat.conjugate()
                
            if bReorientBones:
                (new_bone_size, quat_orient_diff) = calc_bone_rotation(psk_bone, bone_size_choosen, bReorientDirectly, sum_bone_pos)
                # @
                # post_quat = psk_bone.orig_quat.conjugated() * quat_orient_diff

                post_quat = quat_orient_diff
                post_quat.rotate( psk_bone.orig_quat.conjugated() )
            else:
                post_quat = psk_bone.orig_quat.conjugated()
            
            # only length of this vector is matter?
            edit_bone.tail = Vector(( 0.0, new_bone_size, 0.0))

            # @
            # edit_bone.matrix = psk_bone.mat_world * post_quat.to_matrix().to_4x4()

            m = post_quat.copy()
            m.rotate( psk_bone.mat_world )

            m = m.to_matrix().to_4x4()
            m.translation = psk_bone.mat_world.translation

            edit_bone.matrix = m
            
            
            # some dev code...
            #### FINAL
            # post_quat = psk_bone.orig_quat.conjugated() * quat_diff
            # edit_bone.matrix = psk_bone.mat_world * test_quat.to_matrix().to_4x4()
            # edit_bone["post_quat"] = test_quat
            #### 

            # edit_bone["post_quat"] = Quaternion((1,0,0,0))
            # edit_bone.matrix = psk_bone.mat_world* psk_bone.rot
     

            # if edit_bone.parent:
              # edit_bone.matrix = edit_bone.parent.matrix * psk_bone.trans * (psk_bone.orig_quat.conjugated().to_matrix().to_4x4())
              # edit_bone.matrix = edit_bone.parent.matrix * psk_bone.trans * (test_quat.to_matrix().to_4x4())
            # else:
              # edit_bone.matrix = psk_bone.orig_quat.to_matrix().to_4x4()
              
              
            # save bindPose information for .psa import
            # dev
            edit_bone["orig_quat"] = psk_bone.orig_quat
            edit_bone["orig_loc"]  = psk_bone.orig_loc
            edit_bone["post_quat"] = post_quat

            '''
            bone = edit_bone
            if psk_bone.parent is not None:
                orig_loc  =  bone.matrix.translation - bone.parent.matrix.translation
                orig_loc.rotate( bone.parent.matrix.to_quaternion().conjugated() )

                
                orig_quat = bone.matrix.to_quaternion()
                orig_quat.rotate( bone.parent.matrix.to_quaternion().conjugated()  )
                orig_quat.conjugate()

                if orig_quat.dot( psk_bone.orig_quat ) < 0.95:
                    print(bone.name, psk_bone.orig_quat, orig_quat, orig_quat.dot( psk_bone.orig_quat ))
                    print('parent:', bone.parent.matrix.to_quaternion(), bone.parent.matrix.to_quaternion().rotation_difference(bone.matrix.to_quaternion()) )


                if (psk_bone.orig_loc - orig_loc).length > 0.02:
                    print(bone.name, psk_bone.orig_loc, orig_loc, (psk_bone.orig_loc - orig_loc).length)
            '''
    utils_set_mode('OBJECT')
         
    #==================================================================================================
    # Weights
    if bImportmesh: 
    
        vertices_total = len(Vertices)
        
        for ( _, PointIndex, BoneIndex ) in Weights:
            if PointIndex < vertices_total: # can it be not?
                psk_bones[BoneIndex].have_weight_data = True
            # else:
                # print(psk_bones[BoneIndex].name, 'for other mesh',PointIndex ,vertices_total)
                
            #print("weight:", PointIndex, BoneIndex, Weight)
        # Weights.append(None)
        # print(Weights.count(None))
        
        
    # Original vertex colorization code
    '''
    # Weights.sort( key = lambda wgh: wgh[0])
    if bImportmesh:
        VtxCol = []
        bones_count = len(psk_bones)
        for x in range(bones_count):
            #change the overall darkness of each material in a range between 0.1 and 0.9
            tmpVal = ((float(x) + 1.0) / bones_count * 0.7) + 0.1
            tmpVal = int(tmpVal * 256)
            tmpCol = [tmpVal, tmpVal, tmpVal, 0]
            #Change the color of each material slightly
            if x % 3 == 0:
                if tmpCol[0] < 128:
                    tmpCol[0] += 60
                else:
                    tmpCol[0] -= 60
            if x % 3 == 1:
                if tmpCol[1] < 128:
                    tmpCol[1] += 60
                else:
                    tmpCol[1] -= 60
            if x % 3 == 2:
                if tmpCol[2] < 128:
                    tmpCol[2] += 60
                else:
                    tmpCol[2] -= 60
            #Add the material to the mesh
            VtxCol.append(tmpCol)
            
    for x in range(len(Tmsh.faces)):
        for y in range(len(Tmsh.faces[x].v)):
            #find v in Weights[n][0]
            findVal = Tmsh.faces[x].v[y].index
            n = 0
            while findVal != Weights[n][0]:
                n = n + 1
            TmpCol = VtxCol[Weights[n][1]]
            #check if a vertex has more than one influence
            if n != len(Weights) - 1:
                if Weights[n][0] == Weights[n + 1][0]:
                    #if there is more than one influence, use the one with the greater influence
                    #for simplicity only 2 influences are checked, 2nd and 3rd influences are usually very small
                    if Weights[n][2] < Weights[n + 1][2]:
                        TmpCol = VtxCol[Weights[n + 1][1]]
        Tmsh.faces[x].col.append(NMesh.Col(TmpCol[0], TmpCol[1], TmpCol[2], 0))
    '''

    #===================================================================================================
    # UV. Setup.
    
    if bImportmesh:
        # Trick! Create UV maps BEFORE mesh and get (0,0) coordinates for free!
        #   ...otherwise UV coords will be copied from active, or calculated from mesh...
        
        if bSpltiUVdata:
            
            for i in range(len(uv_mat_ids)):
                get_uv_layers(mesh_data).new(name = NAME_UV_PREFIX + str(i))
                
        else:
        
            get_uv_layers(mesh_data).new(name = NAME_UV_PREFIX+"_SINGLE")
        
        
        for counter, uv_data in enumerate(Extrauvs):
            
            if len(mesh_data.uv_layers) < MAX_UVS:
                    
                get_uv_layers(mesh_data).new(name = "EXTRAUVS"+str(counter))
                
            else:
            
                Extrauvs.remove(uv_data)
                print('Extra UV layer %s is ignored. Re-import without "Split UV data".' % counter)
          
    #================================================================================================== 
    # Mesh. Build.
    
        mesh_data.from_pydata(Vertices,[],Faces)

    #==================================================================================================
    # Vertex Normal. Set.

        if Normals is not None:
            mesh_data.polygons.foreach_set("use_smooth", [True] * len(mesh_data.polygons))
            mesh_data.normals_split_custom_set_from_vertices(Normals)
            mesh_data.use_auto_smooth = True
                
    #===================================================================================================
    # UV. Set.
    
    if bImportmesh:

        for face in mesh_data.polygons:
            face.material_index = UV_by_face[face.index][1]

        uv_layers = mesh_data.uv_layers
        
        if not bSpltiUVdata:
           uvLayer = uv_layers[0]
           
        # per face
        # for faceIdx, (faceUVs, faceMatIdx, _, _, wmidx) in enumerate(UV_by_face):
        for faceIdx, (faceUVs, faceMatIdx, WedgeMatIds) in enumerate(UV_by_face):
        
            # per vertex
            for vertN, uv in enumerate(faceUVs):
                loopId = faceIdx * 3 + vertN
                
                if bSpltiUVdata:
                    uvLayer = uv_layers[WedgeMatIds[vertN]]
                    
                uvLayer.data[loopId].uv = uv

    #==================================================================================================
    # VertexColors
    
        if VertexColors is not None:
        
            vtx_color_layer = mesh_data.vertex_colors.new(name = "PSKVTXCOL_0", do_init = False)
            
            pervertex = [None] * len(Vertices)
            
            for counter, (vertexid,_,_,_) in enumerate(Wedges):
            
                # Is it possible ?
                if (pervertex[vertexid] is not None) and (pervertex[vertexid] != VertexColors[counter]):
                    print('Not equal vertex colors. ', vertexid, pervertex[vertexid], VertexColors[counter])
                
                pervertex[vertexid] = VertexColors[counter]
            
            
            for counter, loop in enumerate(mesh_data.loops):
            
                color = pervertex[ loop.vertex_index ]
                
                if color is None:
                    vtx_color_layer.data[ counter ].color = (1.,1.,1.,1.)
                else:
                    if bToSRGB:
                        vtx_color_layer.data[ counter ].color = (
                            color_linear_to_srgb(color[0] / 255),
                            color_linear_to_srgb(color[1] / 255),
                            color_linear_to_srgb(color[2] / 255),
                            color[3] / 255
                        )
                    else:
                        vtx_color_layer.data[ counter ].color = (
                            color[0] / 255,
                            color[1] / 255,
                            color[2] / 255,
                            color[3] / 255
                        )
                        
    #===================================================================================================
    # Extra UVs. Set.
        
        # for counter, uv_data in enumerate(Extrauvs):
        
        #     uvLayer = mesh_data.uv_layers[ counter - len(Extrauvs) ]
            
        #     for uv_index, uv_coords in enumerate(uv_data):
            
        #         uvLayer.data[uv_index].uv = (uv_coords[0], 1.0 - uv_coords[1])


        for counter, uv_data in enumerate(Extrauvs):

            uvLayer = mesh_data.uv_layers[ counter - len(Extrauvs) ]

            for faceIdx, (WedgeIdx3,WedgeIdx2,WedgeIdx1) in enumerate(WedgeIdx_by_faceIdx):
                
                # equal to gltf
                uvLayer.data[faceIdx*3  ].uv = (uv_data[WedgeIdx2][0], 1.0 - uv_data[WedgeIdx2][1])
                uvLayer.data[faceIdx*3+1].uv = (uv_data[WedgeIdx1][0], 1.0 - uv_data[WedgeIdx1][1])
                uvLayer.data[faceIdx*3+2].uv = (uv_data[WedgeIdx3][0], 1.0 - uv_data[WedgeIdx3][1])
                # uvLayer.data[faceIdx*3  ].uv = (uv_data[WedgeIdx3][0], 1.0 - uv_data[WedgeIdx3][1])
                # uvLayer.data[faceIdx*3+1].uv = (uv_data[WedgeIdx2][0], 1.0 - uv_data[WedgeIdx2][1])
                # uvLayer.data[faceIdx*3+2].uv = (uv_data[WedgeIdx1][0], 1.0 - uv_data[WedgeIdx1][1])
        
                
    #===================================================================================================
    # Mesh. Vertex Groups. Bone Weights.
        
        for psk_bone in psk_bones:
            if psk_bone.have_weight_data:
                psk_bone.vertex_group = mesh_obj.vertex_groups.new(name = psk_bone.name)
            # else:
                # print(psk_bone.name, 'have no influence on this mesh')
        
        for weight, vertex_id, bone_index_w in filter(None, Weights):
            psk_bones[bone_index_w].vertex_group.add((vertex_id,), weight, 'ADD')
        
    
    #===================================================================================================
    # Skeleton. Colorize.
    
    if bImportbone:
    
        bone_group_unused = armature_obj.pose.bone_groups.new(name = "Unused bones")
        bone_group_unused.color_set = 'THEME14'

        bone_group_nochild = armature_obj.pose.bone_groups.new(name = "No children")
        bone_group_nochild.color_set = 'THEME03'

        armature_data.show_group_colors = True

        for psk_bone in psk_bones:
        
            pose_bone = armature_obj.pose.bones[psk_bone.name]
            
            if psk_bone.have_weight_data:
            
                if len(psk_bone.children) == 0:
                    pose_bone.bone_group = bone_group_nochild
                    
            else:
                pose_bone.bone_group = bone_group_unused
            
                    
    #===================================================================================================
    # Final
    
    if bImportmesh:
    
        util_obj_link(context, mesh_obj)
        util_select_all(False)
        
            
        if not bImportbone:   
        
            util_obj_select(context, mesh_obj)
            util_obj_set_active(context, mesh_obj)
            
        else:
            # select_all(False)
            util_obj_select(context, armature_obj)
            
            # parenting mesh to armature object
            mesh_obj.parent = armature_obj
            mesh_obj.parent_type = 'OBJECT'
            
            # add armature modifier
            blender_modifier = mesh_obj.modifiers.new( armature_obj.data.name, type = 'ARMATURE')
            blender_modifier.show_expanded = False
            blender_modifier.use_vertex_groups = True
            blender_modifier.use_bone_envelopes = False
            blender_modifier.object = armature_obj
            
            # utils_set_mode('OBJECT')
            # select_all(False)
            util_obj_select(context, armature_obj)
            util_obj_set_active(context, armature_obj)
    
    # print("Done: %f sec." % (time.process_time() - ref_time))
    utils_set_mode('OBJECT')
    return True


class class_psa_bone:
    name = ""
    
    parent = None
    
    fcurve_loc_x = None
    fcurve_loc_y = None
    fcurve_loc_z = None
    fcurve_quat_x = None
    fcurve_quat_y = None
    fcurve_quat_z = None
    fcurve_quat_w = None
    
    post_quat = None
    orig_quat = None
    orig_loc = None


    
def blen_get_armature_from_selection():

  armature_obj = None
  
  for obj in bpy.data.objects:
      if obj.type == 'ARMATURE' and obj_select_get(obj):
          armature_obj = obj
          break
          
  if armature_obj is None:  
    for obj in bpy.data.objects:
      if obj.type == 'MESH' and obj_select_get(obj):
        for modifier in obj.modifiers:
          if modifier.type == 'ARMATURE':
            armature_obj = modifier.object
            break
    
  return armature_obj
  
    
def psaimport(filepath,
        context = None,
        oArmature = None,
        bFilenameAsPrefix = False,
        bActionsToTrack = False,
        first_frames = 0,
        bDontInvertRoot = True,
        bUpdateTimelineRange = False,
        bRotationOnly = False,
        bScaleDown = True,
        fcurve_interpolation = 'LINEAR',
        # error_callback = __pass
        error_callback = print
        ):
    """Import animation data from 'filepath' using 'oArmature'
    
    Args:
        first_frames: (0 - import all)
            Import only 'first_frames' from each action
          
        bActionsToTrack:
            Put all imported actions in one NLAtrack.
          
        oArmature:
            Skeleton used to calculate keyframes
    """
    print ("-----------------------------------------------")
    print ("---------EXECUTING PSA PYTHON IMPORTER---------")
    print ("-----------------------------------------------")
    
    file_ext = 'psa'
    try:
        psafile = open(filepath, 'rb')
    except IOError:
        error_callback('Error while opening file for reading:\n  "'+filepath+'"')
        return False
   
    print ("Importing file: ", filepath)
    
    
    if not context:
        context = bpy.context
        
    armature_obj = oArmature
    
    if armature_obj is None:  
        armature_obj = blen_get_armature_from_selection()
        if armature_obj is None:
            error_callback("No armature selected.")
            return False


    chunk_id = None
    chunk_type = None
    chunk_datasize = None
    chunk_datacount = None
    chunk_data = None

    def read_chunk():
        nonlocal chunk_id, chunk_type,\
                 chunk_datasize, chunk_datacount,\
                 chunk_data

        (chunk_id, chunk_type,
         chunk_datasize, chunk_datacount) = unpack('20s3i', psafile.read(32))
        
        chunk_data = psafile.read(chunk_datacount * chunk_datasize)
    #============================================================================================== 
    # General Header
    #============================================================================================== 
    read_chunk()
    
    if not util_is_header_valid(filepath, file_ext, chunk_id, error_callback):
        return False
    
    #============================================================================================== 
    # Bones (FNamedBoneBinary)
    #============================================================================================== 
    read_chunk()
        
    psa_bones = {}
    
    def new_psa_bone(bone, pose_bone):
        psa_bone = class_psa_bone()
        
        psa_bones[pose_bone.name] = psa_bone
        
        psa_bone.name = pose_bone.name
        
        psa_bone.pose_bone = pose_bone
        
        if bone.parent != None:
            # does needed parent bone was added from psa file
            if bone.parent.name in psa_bones:
                psa_bone.parent = psa_bones[bone.parent.name]
            # no. armature doesnt match
            else:
                psa_bone.parent = None
        # else:
            # psa_bone.parent = None

        # brute fix for non psk skeletons
        if bone.get('orig_quat') is None:

            if bone.parent != None:
                
                psa_bone.orig_loc  =  bone.matrix_local.translation - bone.parent.matrix_local.translation
                psa_bone.orig_loc.rotate( bone.parent.matrix_local.to_quaternion().conjugated() )

                psa_bone.orig_quat = bone.matrix_local.to_quaternion()
                psa_bone.orig_quat.rotate( bone.parent.matrix_local.to_quaternion().conjugated()  )
                psa_bone.orig_quat.conjugate()
            else:
                psa_bone.orig_loc  = bone.matrix_local.translation.copy()
                psa_bone.orig_quat = bone.matrix_local.to_quaternion()

            psa_bone.post_quat = psa_bone.orig_quat.conjugated()
        else:
            psa_bone.orig_quat = Quaternion(bone['orig_quat'])
            psa_bone.orig_loc  =     Vector(bone['orig_loc'])
            psa_bone.post_quat = Quaternion(bone['post_quat'])

        return psa_bone
        
    #Bones Data
    BoneIndex2Name = [None] * chunk_datacount
    BoneNotFoundList = []
    BonesWithoutAnimation = []
    PsaBonesToProcess = [None] * chunk_datacount
    BonePsaImportedNames = []

    # printlog("Name\tFlgs\tNumChld\tPrntIdx\tQx\tQy\tQz\tQw\tLocX\tLocY\tLocZ\tLength\tXSize\tYSize\tZSize\n")

    
    # for case insensetive comparison
    # key = lowered name
    # value = orignal name
    skeleton_bones_lowered = {}
    
    for blender_bone_name in armature_obj.data.bones.keys():
      skeleton_bones_lowered[blender_bone_name.lower()] = blender_bone_name

    
    for counter in range(chunk_datacount):
        
        # tPrntIdx is -1 for parent; and 0 for other; no more useful data
        # indata = unpack_from('64s3i11f', chunk_data, chunk_datasize * counter)
        (indata) = unpack_from('64s56x', chunk_data, chunk_datasize * counter)
        in_name = util_bytes_to_str(indata[0])
        # bonename = util_bytes_to_str(indata[0]).upper()
        
        in_name_lowered = in_name.lower()
        if in_name_lowered in skeleton_bones_lowered:
            orig_name = skeleton_bones_lowered[in_name_lowered]
            
            count_duplicates = BonePsaImportedNames.count( in_name_lowered )

            if count_duplicates > 0:

                duplicate_name_numbered = in_name_lowered + ('.%03d' % count_duplicates)

                # print('Dup:', in_name_lowered, '~',duplicate_name_numbered)

                # Skeleton have duplicate name too?
                if duplicate_name_numbered in skeleton_bones_lowered:
                    orig_name = orig_name + ('.%03d' % count_duplicates)
                else:
                    # Skip animation import for that bone
                    print(" PSK do not have numbered duplicate name(but PSA have!):", duplicate_name_numbered)
                    BonePsaImportedNames.append(in_name_lowered)
                    continue
                    
                
            # use a skeleton bone name 
            BoneIndex2Name[counter] = orig_name
            PsaBonesToProcess[counter] = new_psa_bone(armature_obj.data.bones[orig_name], 
                                                    armature_obj.pose.bones[orig_name])
            BonePsaImportedNames.append(in_name_lowered)
        else:
            # print("Can't find the bone:", orig_name, in_name_lowered)
            BoneNotFoundList.append(counter)
            
    
    if len(psa_bones) == 0:
        error_callback('No bone was match!\nSkip import!')
        return False
    
    # does anyone care?
    for blender_bone_name in armature_obj.data.bones.keys():
        if BoneIndex2Name.count(blender_bone_name) == 0:
            BonesWithoutAnimation.append(blender_bone_name)
            
    if len(BoneNotFoundList) > 0:
        print('PSA have data for more bones: %i.' % len(BoneNotFoundList))
      
    if len(BonesWithoutAnimation) > 0:
        print('PSA do not have data for %i bones:\n' % len(BonesWithoutAnimation), ', '.join(BonesWithoutAnimation))
    #============================================================================================== 
    # Animations (AniminfoBinary)
    #============================================================================================== 
    read_chunk()

    Raw_Key_Nums = 0
    Action_List = [None] * chunk_datacount
    
    for counter in range(chunk_datacount):
        (action_name_raw,        #0
         group_name_raw,         #1
         Totalbones,             #2
         RootInclude,            #3
         KeyCompressionStyle,    #4
         KeyQuotum,              #5
         KeyReduction,           #6
         TrackTime,              #7
         AnimRate,               #8
         StartBone,              #9
         FirstRawFrame,          #10
         NumRawFrames            #11
        ) = unpack_from('64s64s4i3f3i', chunk_data, chunk_datasize * counter)
                       
        action_name = util_bytes_to_str( action_name_raw )
        group_name = util_bytes_to_str( group_name_raw  )

        Raw_Key_Nums += Totalbones * NumRawFrames
        Action_List[counter] = ( action_name, group_name, Totalbones, NumRawFrames)
        
    #============================================================================================== 
    # Raw keys (VQuatAnimKey) 3f vec, 4f quat, 1f time
    #============================================================================================== 
    read_chunk()
    
    if(Raw_Key_Nums != chunk_datacount):
        error_callback(
                'Raw_Key_Nums Inconsistent.'
                '\nData count found: '+chunk_datacount+
                '\nRaw_Key_Nums:' + Raw_Key_Nums
                )
        return False

    Raw_Key_List = [None] * chunk_datacount
    
    unpack_data = Struct('3f4f4x').unpack_from
 
    for counter in range(chunk_datacount):
        pos = Vector()
        quat = Quaternion()
        
        ( pos.x,  pos.y,  pos.z,
         quat.x, quat.y, quat.z, quat.w
        ) = unpack_data( chunk_data, chunk_datasize * counter)
        
        if bScaleDown:
            Raw_Key_List[counter] = (pos * 0.01, quat)
        else:
            Raw_Key_List[counter] = (pos, quat)
    
    psafile.close()
    
    utils_set_mode('OBJECT')

    # index of current frame in raw input data
    raw_key_index = 0
    
    util_obj_set_active(context, armature_obj)
    
    gen_name_part = util_gen_name_part(filepath)
    
    armature_obj.animation_data_create()
    
    if bActionsToTrack:
        nla_track = armature_obj.animation_data.nla_tracks.new()
        nla_track.name = gen_name_part
        nla_stripes = nla_track.strips

        nla_track_last_frame = 0

        if len(armature_obj.animation_data.nla_tracks) > 0:
            for track in armature_obj.animation_data.nla_tracks:
                if len(track.strips) > 0:
                    if track.strips[-1].frame_end > nla_track_last_frame:
                        nla_track_last_frame = track.strips[-1].frame_end

    
    is_first_action = True
    first_action = None
        
    for counter, (Name, Group, Totalbones, NumRawFrames) in enumerate(Action_List):
        ref_time = time.process_time()
        
        if Group != 'None':
            Name = "(%s) %s" % (Group,Name)
        if bFilenameAsPrefix:
            Name = "(%s) %s" % (gen_name_part, Name)
            
        action = bpy.data.actions.new(name = Name)
        
        # force print usefull information to console(due to possible long execution)
        print("Action {0:>3d}/{1:<3d} frames: {2:>4d} {3}".format(
                counter+1, len(Action_List), NumRawFrames, Name)
              )
        
        if first_frames > 0:
            maxframes = first_frames
            keyframes = min(first_frames, NumRawFrames)
            #dev
            # keyframes += 1
        else:
            maxframes = 99999999
            keyframes = NumRawFrames
            
        # create all fcurves(for all bones) for an action
        # for pose_bone in armature_obj.pose.bones:
        for psa_bone in PsaBonesToProcess:
            if psa_bone is None:
                continue
            pose_bone = psa_bone.pose_bone
            
            data_path = pose_bone.path_from_id("rotation_quaternion")
            psa_bone.fcurve_quat_w = action.fcurves.new(data_path, index = 0)
            psa_bone.fcurve_quat_x = action.fcurves.new(data_path, index = 1)
            psa_bone.fcurve_quat_y = action.fcurves.new(data_path, index = 2)
            psa_bone.fcurve_quat_z = action.fcurves.new(data_path, index = 3)
        
            if not bRotationOnly:
                data_path = pose_bone.path_from_id("location")
                psa_bone.fcurve_loc_x = action.fcurves.new(data_path, index = 0)
                psa_bone.fcurve_loc_y = action.fcurves.new(data_path, index = 1)
                psa_bone.fcurve_loc_z = action.fcurves.new(data_path, index = 2)
            
            # 1. Pre-add keyframes! \0/
            # 2. Set data: keyframe_points[].co[0..1]
            # 3. If 2 is not done, do 4: (important!!!)
            # 4. "Apply" data: fcurve.update()
            #      # added keyframes points by default is breaking fcurve somehow
            #      # bcs they are all at the same position?
            psa_bone.fcurve_quat_w.keyframe_points.add(keyframes)
            psa_bone.fcurve_quat_x.keyframe_points.add(keyframes)
            psa_bone.fcurve_quat_y.keyframe_points.add(keyframes)
            psa_bone.fcurve_quat_z.keyframe_points.add(keyframes)

            if not bRotationOnly:
                psa_bone.fcurve_loc_x.keyframe_points.add(keyframes) 
                psa_bone.fcurve_loc_y.keyframe_points.add(keyframes) 
                psa_bone.fcurve_loc_z.keyframe_points.add(keyframes) 
            
        for i in range(0,min(maxframes, NumRawFrames)):
            # raw_key_index+= Totalbones * 5 #55
            for j in range(Totalbones):
                if j in BoneNotFoundList:
                    raw_key_index += 1
                    continue
                
                psa_bone = PsaBonesToProcess[j]
                # pose_bone = psa_bone.pose_bone
                
                p_pos = Raw_Key_List[raw_key_index][0]
                p_quat = Raw_Key_List[raw_key_index][1]
                
                # @
                # if psa_bone.parent:
                    # quat = (p_quat * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)
                # else:
                #     if bDontInvertRoot:
                #         quat = (p_quat.conjugated() * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)
                #     else:
                        # quat = (p_quat * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)

                q = psa_bone.post_quat.copy()
                q.rotate( psa_bone.orig_quat )

                quat = q

                q = psa_bone.post_quat.copy()

                if psa_bone.parent == None and bDontInvertRoot:
                    q.rotate( p_quat.conjugated() )
                else:
                    q.rotate( p_quat )

                quat.rotate( q.conjugated() )
                        
                # @
                # loc = psa_bone.post_quat.conjugated() * p_pos -  psa_bone.post_quat.conjugated() * psa_bone.orig_loc
                
                if not bRotationOnly:
                    loc = (p_pos - psa_bone.orig_loc)
                    # "edit bone" location is in "parent space"
                    # but "pose bone" location is in "local space(bone)"
                    # so we need to transform from parent(edit_bone) to local space (pose_bone)
                    loc.rotate( psa_bone.post_quat.conjugated() )
                    
                # if not bRotationOnly:
                    # loc = (p_pos - psa_bone.orig_loc)
                    # if psa_bone.parent is not None:
                        # q = psa_bone.parent.post_quat.copy()
                        # q.rotate( psa_bone.parent.orig_quat )
                        # print(q)
                        # loc.rotate( psa_bone.parent.post_quat.conjugated() )
                        # loc.rotate( q.conjugated() )
                        # loc.rotate( q )
                        # pass
                
                # quat = p_quat.conjugated()
                # quat = p_quat
                # quat.rotate( psa_bone.orig_quat.conjugated() )
                # quat = Quaternion()
                # loc = -p_pos
                # loc = (p_pos - psa_bone.orig_loc)
                # loc = Vector()
                # loc.rotate( psa_bone.post_quat.conjugated() )

                # Set it?
                # pose_bone.rotation_quaternion = quat
                # pose_bone.location = loc

                    # pose_bone.rotation_quaternion = orig_rot.conjugated()
                    # pose_bone.location = p_pos - (pose_bone.bone.matrix_local.translation - pose_bone.bone.parent.matrix_local.translation)
                
                ##### Works + post_quat (without location works)
                # quat = (p_quat * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)
                # loc = psa_bone.post_quat.conjugated() * (p_pos - psa_bone.orig_loc)

                
                psa_bone.fcurve_quat_w.keyframe_points[i].co = i, quat.w
                psa_bone.fcurve_quat_x.keyframe_points[i].co = i, quat.x
                psa_bone.fcurve_quat_y.keyframe_points[i].co = i, quat.y
                psa_bone.fcurve_quat_z.keyframe_points[i].co = i, quat.z
                
                psa_bone.fcurve_quat_w.keyframe_points[i].interpolation = fcurve_interpolation
                psa_bone.fcurve_quat_x.keyframe_points[i].interpolation = fcurve_interpolation
                psa_bone.fcurve_quat_y.keyframe_points[i].interpolation = fcurve_interpolation
                psa_bone.fcurve_quat_z.keyframe_points[i].interpolation = fcurve_interpolation
                
                
                if not bRotationOnly:
                    psa_bone.fcurve_loc_x.keyframe_points[i].co = i, loc.x
                    psa_bone.fcurve_loc_y.keyframe_points[i].co = i, loc.y
                    psa_bone.fcurve_loc_z.keyframe_points[i].co = i, loc.z
                    
                    psa_bone.fcurve_loc_x.keyframe_points[i].interpolation = fcurve_interpolation
                    psa_bone.fcurve_loc_y.keyframe_points[i].interpolation = fcurve_interpolation
                    psa_bone.fcurve_loc_z.keyframe_points[i].interpolation = fcurve_interpolation
                
                # Old path. Slower.
                # psa_bone.fcurve_quat_w.keyframe_points.insert(i,quat.w,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # psa_bone.fcurve_quat_x.keyframe_points.insert(i,quat.x,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # psa_bone.fcurve_quat_y.keyframe_points.insert(i,quat.y,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # psa_bone.fcurve_quat_z.keyframe_points.insert(i,quat.z,{'NEEDED','FAST'}).interpolation = fcurve_interpolation

                # psa_bone.fcurve_loc_x.keyframe_points.insert(i,loc.x,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # psa_bone.fcurve_loc_y.keyframe_points.insert(i,loc.y,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # psa_bone.fcurve_loc_z.keyframe_points.insert(i,loc.z,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                raw_key_index += 1
            
            # on first frame
            # break
        raw_key_index += (NumRawFrames-min(maxframes,NumRawFrames)) * Totalbones

        # Add action to tail of the nla track
        if bActionsToTrack:
            
            if len(nla_track.strips) == 0:
                strip = nla_stripes.new(Name, nla_track_last_frame, action)
            else:
                strip = nla_stripes.new(Name, nla_stripes[-1].frame_end, action)

            # Do not pollute track. Makes other tracks 'visible' through 'empty space'.
            strip.extrapolation = 'NOTHING'

            nla_track_last_frame += NumRawFrames

        if is_first_action:
            first_action = action
            is_first_action = False
            
        print("Done: %f sec." % (time.process_time() - ref_time))
        # break on first animation set
        # break
        
    scene = util_get_scene(context)

    if not bActionsToTrack:
        if not scene.is_nla_tweakmode:
            armature_obj.animation_data.action = first_action
    
    if bUpdateTimelineRange:

        scene.frame_start = 0

        if bActionsToTrack:
            scene.frame_end = sum(frames for _, _, _, frames in Action_List) - 1
        else:
            scene.frame_end = max(frames for _, _, _, frames in Action_List) - 1


    util_select_all(False)
    util_obj_select(context, armature_obj)
    util_obj_set_active(context, armature_obj)
    
    # 2.8 crashes
    # scene.frame_set(0)

class PSKPSA_OT_show_message(bpy.types.Operator):
    bl_idname = "pskpsa.message"
    bl_label = "PSA/PSK"
    bl_options = {'REGISTER', 'INTERNAL'}

    message : StringProperty(default = 'Message')

    lines = []
    line0 = None
    def execute(self, context):
        self.lines = self.message.split("\n")
        maxlen = 0
        for line in self.lines:
            if len(line) > maxlen:
                maxlen = len(line)
                
        print(self.message)
            
        self.report({'WARNING'}, self.message)
        return {'FINISHED'}
        
    def invoke(self, context, event):
        self.lines = self.message.split("\n")
        maxlen = 0
        for line in self.lines:
            if len(line) > maxlen:
                maxlen = len(line)
                
        self.line0 = self.lines.pop(0)
        
        return context.window_manager.invoke_props_dialog(self, width = 100 + 6*maxlen)
      
    def cancel(self, context):
        # print('cancel')
        self.execute(self)
        
    def draw(self, context):
        layout = self.layout
        sub = layout.column()
        sub.label(text = self.line0, icon = 'ERROR')

        for line in self.lines:
            sub.label(text = line)
        

#properties for panels, and Operator.
class ImportProps():

    fBonesize : FloatProperty(
            name = "Alt. bone length.",
            description = "Bone length will be set to this value IF it's less than [Corrected avg. bone length]\nBone length = min( <this value> , [Corrected avg. bone length] )",
            default = 5.0, min = 0.01, max = 50, step = 0.3, precision = 2,
            )
    fBonesizeRatio : FloatProperty(
            name = "Bone length ratio",
            description = "Bone length will be set to this value IF it's less than  Corrected avg. bone length = <this value> * [calculated average bone length]",
            default = 0.4, min = 0.1, max = 4, step = 0.05, precision = 2,
            )
    bSpltiUVdata : BoolProperty(
            name = "Split UV data",
            description = "Try to place UV points(coords) to different UV maps, according to material index of the Wedge(UV-Vertex-MateralIndex)."\
                    "\n * Each UVmap will still have the same amount of points(data/coords)."\
                    "\n * Blender can have only 8 UVs per mesh. So it is not enough to have different UVmap for each material.",
            default = False,
            )
    bReorientBones : BoolProperty(
            name = "Reorient bones",
            description = "Bones will be axis-aligned to children.",
            default = False,
            )
    bReorientDirectly : BoolProperty(
            name = "Reorient directly",
            description = "Directly to children.\n * Axes will not be preserved.\n * Orphan bones - in direction from parent head\n * With only one non-orphan bone - to that one.",
            default = False,
            ) 
    import_mode : EnumProperty(
            name = "Import mode.",
            items = (('All','All','Import mesh and skeleton'),
                    ('Mesh','Mesh','Import only mesh'),
                    ('Skel','Skel','Import only skeleton'))
            )
    bDontInvertRoot : BoolProperty(
            name = "Don't invert root bone",
            description = " * Used by PSK and PSA.\n * Uncheck it, if skeleton is badly oriented.",
            default = True,
            )
    bFilenameAsPrefix :  BoolProperty(
            name = "Prefix action names",
            description = "Use the filename as a prefix for the action name.",
            default = False,
            )
    bActionsToTrack : BoolProperty(
            name = "All actions to NLA track",
            description = "Add all imported actions to new NLAtrack. One by one.\nStarting at the very end.\nLook at \"Nonlinear Animation\" editor.",
            default = False,
            )
    bUpdateTimelineRange : BoolProperty(
            name = "Update timeline range",
            description = "Set timeline range to match imported action[s] length.\n * If \"All actions to NLA track\" is disabled, range will be set to hold longest action.",
            default = False,
            )
    bRotationOnly : BoolProperty(
            name = "Rotation only",
            description = "Create only rotation keyframes.",
            default = False,
            )
    bScaleDown : BoolProperty(
            name = "Scale down",
            description = " * Used by PSK and PSA.\n * Multiply coordinates by 0.01\n * From \"cm.\" to \"m.\"",
            default = True,
            )
    bToSRGB : BoolProperty(
            name = "sRGB vertex color",
            description = "Apply 'linear RGB -> sRGB' conversion over vertex colors",
            default = True,
            )
            
    def draw_psk(self, context):
        props = bpy.context.scene.pskpsa_import
        layout = self.layout
        layout.prop(props, 'import_mode', expand = True)
        layout.prop(props, 'bReorientBones')
        
        sub = layout.row()
        sub.prop(props, 'bReorientDirectly')
        sub.enabled = props.bReorientBones
        
        # layout.prop(props, 'bDontInvertRoot')
        layout.prop(props, 'bSpltiUVdata')
        sub = layout.row()
        # layout.prop(props, 'bDontInvertRoot', icon = 'ERROR' if props.bDontInvertRoot else 'NONE')
        sub.prop(props, 'bDontInvertRoot')
        if not props.bDontInvertRoot:
            sub.label(text = "", icon = 'ERROR')
            
        layout.prop(props, 'bScaleDown')
        layout.prop(props, 'bToSRGB')
        layout.prop(props, 'fBonesizeRatio')
        layout.prop(props, 'fBonesize')
        
    def draw_psa(self, context):
        props = context.scene.pskpsa_import
        layout = self.layout
        layout.prop(props,'bActionsToTrack')
        layout.prop(props,'bFilenameAsPrefix')
        layout.prop(props,'bUpdateTimelineRange')
        layout.prop(props,'bRotationOnly')
        # layout.prop(props, 'bDontInvertRoot')
        # layout.separator()
   
class PskImportOptions(bpy.types.PropertyGroup, ImportProps):
    pass

def blen_hide_unused(armature_obj, mesh_obj):
    def is_bone_useless(psa_bone):
        is_useless = True
        if len(psa_bone.children) == 0:
            if mesh_obj.vertex_groups.get(psa_bone.name) != None:
                is_useless = False
        else:
            for psa_bone_child in psa_bone.children:
                is_useless = is_bone_useless(psa_bone_child)
                if not is_useless:
                    is_useless = False
                    # break
        # print(psa_bone.name, is_useless)
        if is_useless:
            psa_bone.hide = True
        return is_useless

    # print(armature_obj.data.bones[0].name)
    is_bone_useless(armature_obj.data.bones[0])


class PSKPSA_OT_hide_unused_bones(bpy.types.Operator):
    """Hide bones with no weights and no children\n* Select mesh with armature modifier.\n* ALT + H to reveal (in Pose Mode(CTRL + TAB))"""
    bl_idname = "armature.hide_unused"
    bl_label = "Hide useless bones"
    bl_options = {'UNDO'}

    def execute(self, context):
        mesh_obj = None
        if context.object.type == 'MESH':
            for mod in context.object.modifiers:
                if mod.type == 'ARMATURE' and mod.object:
                    blen_hide_unused(
                        context.object.modifiers[0].object, context.object)
                    return {'FINISHED'}
                else:
                    util_ui_show_msg(
                        "Can't find any situable Armature modifier for selected mesh.")
                        
        elif context.object.type == 'ARMATURE':
            for obj in bpy.context.selected_objects:
                if obj.type == 'MESH':
                  for modifier in obj.modifiers:
                    if modifier.type == 'ARMATURE':
                      if modifier.object == context.object:
                          blen_hide_unused(
                                      context.object, context.object)
                      return {'FINISHED'}
            
        return {'FINISHED'}

        
class IMPORT_OT_psk(bpy.types.Operator, ImportProps):
    
    bl_idname = "import_scene.psk"
    bl_label = "Import PSK"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_options = {'UNDO'}

    filepath : StringProperty(
            subtype = 'FILE_PATH',
            )
    filter_glob : StringProperty(
            default = "*.psk;*.pskx",
            options = {'HIDDEN'},
            )
    files : bpy.props.CollectionProperty(type=bpy.types.OperatorFileListElement, options={'HIDDEN', 'SKIP_SAVE'})
    directory : bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
    
    def draw(self, context):
        self.draw_psk(context)
        # self.layout.prop(context.scene.pskpsa_import, 'bDontInvertRoot')
    # draw = ImportProps.draw_psk
    
    def execute(self, context):
        if not self.filepath:
            raise Exception("filepath not set")
            
        no_errors = True
        
        if not self.directory:
            # possibly excuting from script, 
            # bcs blender will set this value, even for a single file
            
            keywords = self.as_keywords(
                ignore=(
                    "import_mode",
                    "filter_glob",
                    "bFilenameAsPrefix",
                    "bActionsToTrack",
                    "bUpdateTimelineRange",
                    "bRotationOnly",
                    "files", 
                    "directory"
                    )
                )
                
            if self.import_mode == 'Mesh':
                bImportmesh = True
                bImportbone = False
            elif self.import_mode == 'Skel':
                bImportmesh = False
                bImportbone = True
            else:
                bImportmesh = True
                bImportbone = True
            
            # ugly workaround
            keywords["bImportbone"] = bImportbone
            keywords["bImportmesh"] = bImportmesh
            
            no_errors = pskimport( **keywords )
            
        else:        
            props = bpy.context.scene.pskpsa_import
            if props.import_mode == 'Mesh':
                bImportmesh = True
                bImportbone = False
            elif props.import_mode == 'Skel':
                bImportmesh = False
                bImportbone = True
            else:
                bImportmesh = True
                bImportbone = True
            
            
            for _, fileListElement in enumerate(self.files):
                fpath = self.directory + fileListElement.name
                
                no_errors = no_errors and pskimport( 
                            fpath,
                            context = context,
                            bImportmesh = bImportmesh, bImportbone = bImportbone,
                            fBonesize = props.fBonesize,
                            fBonesizeRatio = props.fBonesizeRatio,
                            bSpltiUVdata = props.bSpltiUVdata,
                            bReorientBones = props.bReorientBones,
                            bReorientDirectly = props.bReorientDirectly,
                            bDontInvertRoot = props.bDontInvertRoot,
                            bScaleDown = props.bScaleDown,
                            bToSRGB = props.bToSRGB,
                            error_callback = util_ui_show_msg
                            )

        if not no_errors:
            return {'CANCELLED'}
        else:
            return {'FINISHED'}
    
    def invoke(self, context, event):
        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}

        
class IMPORT_OT_psa(bpy.types.Operator, ImportProps):
    '''Load a skeleton animation from .psa\n * Selected armature will be used.'''
    bl_idname = "import_scene.psa"
    bl_label = "Import PSA"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_options = {'UNDO'}

    filepath : StringProperty(
            subtype = 'FILE_PATH',
            )
    filter_glob : StringProperty(
            default = "*.psa",
            options = {'HIDDEN'},
            )
    files : bpy.props.CollectionProperty(type=bpy.types.OperatorFileListElement, options={'HIDDEN', 'SKIP_SAVE'})
    directory : bpy.props.StringProperty(subtype='FILE_PATH', options={'HIDDEN', 'SKIP_SAVE'})
            
    def draw(self, context):
        self.draw_psa(context)
        self.layout.prop(context.scene.pskpsa_import, 'bDontInvertRoot')
      
    def execute(self, context):
        props = context.scene.pskpsa_import
        
        if not self.directory:
            # possibly excuting from script, 
            # bcs blender will set this value, even for a single file
            psaimport( **(self.as_keywords(
                ignore=(
                    "import_mode",
                    "fBonesize",
                    "fBonesizeRatio",
                    "bSpltiUVdata",
                    "bReorientBones",
                    "bReorientDirectly",
                    "bToSRGB",
                    "filter_glob",
                    "files", 
                    "directory"
                    )
                )) )
            return {'FINISHED'}
        
        for _, fileListElement in enumerate(self.files):
            fpath = self.directory + fileListElement.name
            psaimport(
                fpath,
                context = context,
                bFilenameAsPrefix = props.bFilenameAsPrefix, 
                bActionsToTrack = props.bActionsToTrack, 
                oArmature = blen_get_armature_from_selection(),
                bDontInvertRoot = props.bDontInvertRoot,
                bUpdateTimelineRange = props.bUpdateTimelineRange,
                bRotationOnly = props.bRotationOnly,
                bScaleDown = props.bScaleDown,
                error_callback = util_ui_show_msg
                )
        return {'FINISHED'}
    
    def invoke(self, context, event):
        if blen_get_armature_from_selection() is None:
            util_ui_show_msg('Select an armature.')
            return {'FINISHED'}
        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}

class PSKPSA_PT_import_panel(bpy.types.Panel, ImportProps):
    bl_label = "PSK/PSA Import"
    bl_idname = "VIEW3D_PT_udk_import_280"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "PSK / PSA"
    
    # @classmethod
    # def poll(cls, context):
        # print(context.scene.get('pskpsa_import'),'poll')
        # context.scene.update_tag()
        # context.scene.update()
        # return context.scene.get('pskpsa_import') is not None
          
    def draw(self, context):
        props = context.scene.pskpsa_import
        if props is None:
            self.layout.label(text = "??")
            return
        # return
        layout = self.layout
       
        # layout.label(text = "Mesh and skeleton:")
        layout.operator(IMPORT_OT_psk.bl_idname, icon = 'MESH_DATA')
        self.draw_psk(context)
        # layout.prop(props, 'import_mode',expand = True)
        
        sub = layout.row()
        sub.operator(PSKPSA_OT_hide_unused_bones.bl_idname, icon = 'BONE_DATA')
        sub.enabled = (context.object is not
                       None) and (context.object.type == 'MESH' or context.object.type == 'ARMATURE')
                       
        layout.separator()
        layout.separator()
        # layout.label(text = "Animation:", icon = 'ANIM')
        layout.operator(IMPORT_OT_psa.bl_idname, icon = 'ANIM')
        self.draw_psa(context)

    
def menu_import_draw(self, context):
    self.layout.operator(IMPORT_OT_psk.bl_idname, text = "Skeleton Mesh (.psk)")
    self.layout.operator(IMPORT_OT_psa.bl_idname, text = "Skeleton Anim (.psa)")

classes = (
        IMPORT_OT_psk,
        IMPORT_OT_psa,
        PskImportOptions,
        PSKPSA_PT_import_panel,
        PSKPSA_OT_show_message,
        PSKPSA_OT_hide_unused_bones
    )
    
    
def register():
    from bpy.utils import register_class
    for cls in classes:
        register_class(cls)
        
    bpy.types.TOPBAR_MT_file_import.append(menu_import_draw)

    bpy.types.Scene.pskpsa_import = PointerProperty(type = PskImportOptions)
    
def unregister():
    from bpy.utils import unregister_class
    for cls in classes:
        unregister_class(cls)
        
    bpy.types.TOPBAR_MT_file_import.remove(menu_import_draw)

    del bpy.types.Scene.pskpsa_import
    
if __name__ == "__main__":
    register()

if __name__ == "io_import_scene_unreal_psa_psk_270_dev":
    import pskpsadev
