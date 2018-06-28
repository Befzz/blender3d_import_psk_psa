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
    "version": (2, 7, 2),
    "blender": (2, 64, 0),
    "location": "File > Import > Skeleton Mesh (.psk)/Animation Set (.psa) OR View3D > Tool Shelf (key T) > Misc. tab",
    "description": "Import Skeleton Mesh / Animation Data",
    "warning": "",
    "wiki_url": "http://wiki.blender.org/index.php/Extensions:2.5/Py/"
                "Scripts/Import-Export/Unreal_psk_psa",
    "category": "Import-Export",
    "tracker_url": "https://github.com/Befzz/blender3d_import_psk_psa"
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

# import struct
import time

#DEV
# from mathutils import *
# from math import *

### Cross-version proxy functions 2.8 - 2.7 WiP
# Don't using bpy.context in 2.8
# bpy.app.version: (2,80,17)
v = bpy.app.version
is_blen_280 = (v[0]*1000000 + v[1]*1000 + v[2]) >= 2080017

if is_blen_280:
    def util_obj_link(context, obj):
        # return bpy.context.scene_collection.objects.link(obj)
        # bpy.context.view_layer.collections[0].collection.objects.link(obj)
        # return bpy.context.collection.objects.link(obj)
        # bpy.data.scenes[0].collection.objects.link(obj)
        context.collection.objects.link(obj)

    def util_obj_select(context, obj, action = 'SELECT'):
          # if obj.name in bpy.data.scenes[0].view_layers[0].objects:
          if obj.name in context.view_layer.objects:
              return obj.select_set(action)

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
else:
    def util_obj_link(context, obj):
        bpy.context.scene.objects.link(obj)

    def util_obj_select(context, obj, action = 'SELECT'):
        obj.select = (action == 'SELECT')

    def util_obj_set_active(context, obj):
        bpy.context.scene.objects.active = obj

    def get_uv_layers(mesh_obj):
        return mesh_obj.uv_textures
      
    def util_get_scene(context):
        return bpy.context.scene
        
    def obj_select_get(obj):
        return obj.select
del v

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
            axis_vec = psk_bone.orig_quat * psk_bone.orig_loc
        else:
            # bone Head near parent Head?
            if psk_bone.orig_loc.length < 0.1 * avg_bone_len:
                vec_to_axis_vec(psk_bone.orig_quat.conjugated() * psk_bone.parent.axis_vec, axis_vec)
                # reorient bone to other axis bychanging our base Y vec...
                # this is not tested well
                vecy = Vector((1.0, 0.0, 0.0))
            else:
                vec_to_axis_vec(psk_bone.orig_quat.conjugated() * psk_bone.parent.axis_vec, axis_vec)

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
        
    
def pskimport(filepath,
        context = bpy.context,
        bImportmesh = True,
        bImportbone = True,
        bSpltiUVdata = False,
        fBonesize = 2.0,
        fBonesizeRatio = 0.6,
        bDontInvertRoot = False,
        bReorientBones = False,
        bReorientDirectly = False,
        error_callback = None):
    '''
    Import mesh and skeleton from .psk/.pskx files
    
    Args:
        bReorientBones:
            Axis based bone orientation to children
            
        error_callback:
            Called when importing is failed.
            
            __name__('?', error_callback = lambda msg: print('reason:',msg)
            
    '''
    if not hasattr( error_callback, '__call__'):
        error_callback = __pass
        
    # ref_time = time.process_time()
    if not bImportbone and not bImportmesh:
        error_callback("Nothing to do.\nSet something for import.")
        return False
    file_ext = 'psk'
    
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
    Extrauvs = []
        
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
        nonlocal Faces, UV_by_face
        
        UV_by_face = [None] * chunk_datacount
        Faces = [None] * chunk_datacount
        
        if len(Vertices) > 65536:
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
            Faces[counter] = (vertid2,  vertid1, vertid0)
            
            uv = ( ( u2, 1.0 - v2 ), ( u1, 1.0 - v1 ), ( u0, 1.0 - v0 ) )
            
            # Mapping: FaceIndex <=> UV data <=> FaceMatIndex
            UV_by_face[counter] = (uv, MatIndex, (matid2, matid1, matid0))
            
            
    #==================================================================================================
    # Vertices X | Y | Z
    def read_vertices():
        nonlocal bImportbone, bImportmesh
        
        if not bImportmesh:
            return True
            
        nonlocal Vertices
        
        Vertices = [None] * chunk_datacount
        
        unpack_data = Struct('3f').unpack_from
        
        for counter in range( chunk_datacount ):
            (vec_x, vec_y, vec_z) = unpack_data(chunk_data, counter * chunk_datasize)
            Vertices[counter]  = (vec_x, vec_y, vec_z)
            
            
    #================================================================================================== 
    # Wedges (UV)   VertexId |  U |  V | MatIdx 
    def read_wedges():
    
        nonlocal bImportbone, bImportmesh
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
            unpack_data = Struct('64s3i11f').unpack_from
        else:
            unpack_data = Struct('64s56x').unpack_from
            
        Bones = [None] * chunk_datacount
        
        for counter in range( chunk_datacount ):
            Bones[counter] = unpack_data( chunk_data, chunk_datasize * counter)
          
    
    #================================================================================================== 
    # Influences (Bone Weight) (VRawBoneInfluence) ( Weight | PntIdx | BoneIdx)
    def read_weights():
        # nonlocal Weights, bImportmesh
        nonlocal Weights
        
        if not bImportmesh:
            return True
            
        Weights = [None] * chunk_datacount
        
        unpack_data = Struct('fii').unpack_from
        
        for counter in range(chunk_datacount):
            Weights[counter] = unpack_data(chunk_data, chunk_datasize * counter)
             
    
    #================================================================================================== 
    # Extra UV. U | V
    def read_extrauvs():
        unpack_data = Struct("=2f").unpack_from
        
        uvdata = [None] * chunk_datacount
        
        for counter in range( chunk_datacount ):
            uvdata[counter] = unpack_data(chunk_data, chunk_datasize * counter) 
            
        Extrauvs.append(uvdata)
 
             
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
        'EXTRAUVS': read_extrauvs
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
             vec_x, vec_y, vec_z,                       #8 9 10
             joint_length,                              #11
             scale_x, scale_y, scale_z) in enumerate(Bones):
        
            psk_bone = init_psk_bone(counter, psk_bones, name_raw)
            
            psk_bone.bone_index = counter
            psk_bone.parent_index = ParentIndex
            
            # psk_bone.scale = (scale_x, scale_y, scale_z)

            # store bind pose to make it available for psa-import via CustomProperty of the Blender bone
            psk_bone.orig_quat = Quaternion((quat_w, quat_x, quat_y, quat_z))
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
            psk_bone.mat_world = parent.mat_world_rot.to_4x4()
            psk_bone.mat_world.translation = parent.mat_world.translation + parent.mat_world_rot * psk_bone.orig_loc
            psk_bone.mat_world_rot = parent.mat_world_rot * psk_bone.orig_quat.conjugated().to_matrix()
            
            # psk_bone.mat_world =  ( parent.mat_world_rot.to_4x4() * psk_bone.trans)
            # psk_bone.mat_world.translation += parent.mat_world.translation
            # psk_bone.mat_world_rot = parent.mat_world_rot * psk_bone.orig_quat.conjugated().to_matrix()
            
    
    #==================================================================================================
    # Skeleton. Prepare.
            
        armature_data = bpy.data.armatures.new(gen_names['armature_data'])
        armature_obj = bpy.data.objects.new(gen_names['armature_object'], armature_data)
        # TODO: options for axes and x_ray?
        armature_data.show_axes = False
        armature_data.draw_type = 'STICK'
        armature_obj.show_x_ray = True

        util_obj_link(context, armature_obj)

        util_select_all(False)
        util_obj_select(context, armature_obj)
        util_obj_set_active(context, armature_obj)
        
        utils_set_mode('EDIT')
        
        
        sum_bone_pos /= len(Bones) # average
        sum_bone_pos *= fBonesizeRatio # corrected
        
        bone_size_choosen = max(0.01, min(sum_bone_pos, fBonesize))
        
    #==================================================================================================
    # Skeleton. Build.
        for psk_bone in psk_bones:
            edit_bone = armature_obj.data.edit_bones.new(psk_bone.name)

            armature_obj.data.edit_bones.active = edit_bone

            if psk_bone.parent is not None:
                edit_bone.parent = armature_obj.data.edit_bones[psk_bone.parent.name]
            else:
                if bDontInvertRoot:
                    psk_bone.orig_quat.conjugate()
                
            if bReorientBones:
                (new_bone_size, quat_orient_diff) = calc_bone_rotation(psk_bone, bone_size_choosen, bReorientDirectly, sum_bone_pos)
                post_quat = psk_bone.orig_quat.conjugated() * quat_orient_diff
            else:
                new_bone_size = bone_size_choosen
                post_quat = psk_bone.orig_quat.conjugated()
            
            # only length of this vector is matter?
            edit_bone.tail = Vector(( 0.0, new_bone_size, 0.0))
            # edit_bone.tail = Vector((0.0, 0.0, new_bone_size))
            
            # edit_bone.matrix = psk_bone.mat_world * quat_diff.to_matrix().to_4x4()
            edit_bone.matrix = psk_bone.mat_world * post_quat.to_matrix().to_4x4()
            
            
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
            edit_bone["orig_quat"] = psk_bone.orig_quat
            edit_bone["orig_loc"]  = psk_bone.orig_loc
            edit_bone["post_quat"] = post_quat
            
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

    #===================================================================================================
    # Extra UVs. Set.
        
        for counter, uv_data in enumerate(Extrauvs):
        
            uvLayer = mesh_data.uv_layers[ counter - len(Extrauvs) ]
            
            for uv_index, uv_coords in enumerate(uv_data):
            
                uvLayer.data[uv_index].uv = (uv_coords[0], 1.0 - uv_coords[1])
                
    #===================================================================================================
    # Mesh. Vertex Groups. Bone Weights.
        
        for psk_bone in psk_bones:
            if psk_bone.have_weight_data:
                psk_bone.vertex_group = mesh_obj.vertex_groups.new(psk_bone.name)
            # else:
                # print(psk_bone.name, 'have no influence on this mesh')
        
        for weight, vertex_id, bone_index_w in filter(None, Weights):
            psk_bones[bone_index_w].vertex_group.add((vertex_id,), weight, 'ADD')
        
    
    #===================================================================================================
    # Skeleton. Colorize.
    
    if bImportbone:
    
        bone_group_unused = armature_obj.pose.bone_groups.new("Unused bones")
        bone_group_unused.color_set = 'THEME14'

        bone_group_nochild = armature_obj.pose.bone_groups.new("No children")
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
        context = bpy.context,
        oArmature = None,
        bFilenameAsPrefix = False,
        bActionsToTrack = False,  
        first_frames = 0, 
        bDontInvertRoot = False, 
        fcurve_interpolation = 'LINEAR', 
        error_callback = __pass
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
            
        psa_bone.orig_quat = Quaternion(bone['orig_quat'])
        psa_bone.orig_loc  =     Vector(bone['orig_loc'])
        psa_bone.post_quat = Quaternion(bone['post_quat'])
        return psa_bone
        
    #Bones Data
    BoneIndex2Name = [None] * chunk_datacount
    BoneNotFoundList = []
    BonesWithoutAnimation = []
    PsaBonesToProcess = [None] * chunk_datacount

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
            
            # use a skeleton bone name 
            BoneIndex2Name[counter] = orig_name
            PsaBonesToProcess[counter] = new_psa_bone(armature_obj.data.bones[orig_name], 
                                                      armature_obj.pose.bones[orig_name])
        else:
            # print("Can't find the bone:", bonename)
            BoneNotFoundList.append(counter)
            
    
    if len(psa_bones) == 0:
        error_callback('No bone was match!\nSkip import!')
        return False
    
    # does anyone care?
    for blender_bone_name in armature_obj.data.bones.keys():
        if BoneIndex2Name.count(blender_bone_name) == 0:
            BonesWithoutAnimation.append(blender_bone_name)
            
    if len(BoneNotFoundList) > 0:
      print('Not found bones: %i.' % len(BoneNotFoundList));
      
    if len(BonesWithoutAnimation) > 0:
      print('Bones(%i) without animation data:\n' % len(BonesWithoutAnimation), ', '.join(BonesWithoutAnimation))
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
    else:
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
        
            data_path = pose_bone.path_from_id("location")
            psa_bone.fcurve_loc_x = action.fcurves.new(data_path, index = 0)
            psa_bone.fcurve_loc_y = action.fcurves.new(data_path, index = 1)
            psa_bone.fcurve_loc_z = action.fcurves.new(data_path, index = 2)
            
            # 1. Pre-add keyframes! \0/
            # 2. Set data: keyframe_points[].co[0..1]
            # 3. If 2 is not done, do 4: (important!!!)
            # 4. "Apply" data: fcurve.update()
            # added keyframes points by default is breaking fcurve somehow
            # bcs they are all at the same position?
            psa_bone.fcurve_quat_w.keyframe_points.add(keyframes)
            psa_bone.fcurve_quat_x.keyframe_points.add(keyframes)
            psa_bone.fcurve_quat_y.keyframe_points.add(keyframes)
            psa_bone.fcurve_quat_z.keyframe_points.add(keyframes)

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
                pose_bone = psa_bone.pose_bone
                
                p_pos = Raw_Key_List[raw_key_index][0]
                p_quat = Raw_Key_List[raw_key_index][1]
                
                ##### Worked with no bone rotation
                # quat = p_quat.conjugated() * psa_bone.orig_quat
                # loc = p_pos - psa_bone.orig_loc
                #####
                    

                if psa_bone.parent:
                    ##### Correct
                    # orig_prot = pose_bone.bone.parent.matrix_local.to_3x3().to_quaternion()
                    # orig_rot = pose_bone.bone.matrix_local.to_3x3().to_quaternion()
                    # orig_rot = (orig_prot.conjugated() * orig_rot)
                    ######

                    #### FINAL
                    quat = (p_quat * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)
                    # loc = psa_bone.post_quat.conjugated() * p_pos -  psa_bone.post_quat.conjugated() * psa_bone.orig_loc
                    ####
                else:
                    if bDontInvertRoot:
                        quat = (p_quat.conjugated() * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)
                    else:
                        quat = (p_quat * psa_bone.post_quat).conjugated() * (psa_bone.orig_quat * psa_bone.post_quat)
                        
                loc = psa_bone.post_quat.conjugated() * p_pos -  psa_bone.post_quat.conjugated() * psa_bone.orig_loc
                    
                pose_bone.rotation_quaternion = quat
                pose_bone.location = loc
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
            if nla_track_last_frame == 0:
                nla_stripes.new(Name, 0, action)
            else:
                nla_stripes.new(Name, nla_stripes[-1].frame_end, action)

            nla_track_last_frame += NumRawFrames
        elif is_first_action:
            first_action = action
            is_first_action = False
            
        print("Done: %f sec." % (time.process_time() - ref_time))
        # break on first animation set
        # break
        
    scene = util_get_scene(context)
    
    if not bActionsToTrack:
        if not scene.is_nla_tweakmode:
            armature_obj.animation_data.action = first_action
            
    util_select_all(False)
    util_obj_select(context, armature_obj)
    util_obj_set_active(context, armature_obj)
    
    # 2.8 crashes
    if not is_blen_280:
        scene.frame_set(0)

 
class PSKPSA_OT_show_message(bpy.types.Operator):
    bl_idname = "pskpsa.message"
    bl_label = "PSA/PSK"
    bl_options = {'REGISTER', 'INTERNAL'}

    message = StringProperty(default = 'Message')
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
        sub.label(self.line0, icon = 'ERROR')

        for line in self.lines:
            sub.label(line)
    
#properties for panels, and Operator.
class ImportProps():
    fBonesize = FloatProperty(
            name = "Orphan bone length",
            description = "Maximum orphan bone length",
            default = 5.0, min = 0.1, max = 50, step = 0.3, precision = 2,
            )
    fBonesizeRatio = FloatProperty(
            name = "Bone length ratio",
            description = "Bone length = [average bone length] * [this value]",
            default = 0.6, min = 0.1, max = 4, step = 0.05, precision = 2,
            )
    bSpltiUVdata = BoolProperty(
            name = "Split UV data",
            description = "Try to place UV points(coords) to different UV maps, according to material index of the Wedge(UV-Vertex-MateralIndex)."\
                    "\n * Each UVmap will still have the same amount of points(data/coords)."\
                    "\n * Blender can have only 8 UVs per mesh. So it is not enough to have different UVmap for each material.",
            default = False,
            )
    bReorientBones = BoolProperty(
            name = "Reorient bones",
            description = "Bones will be axis-aligned to children.",
            default = False,
            )
    bReorientDirectly = BoolProperty(
            name = "Reorient directly",
            description = "Directly to children.\n * Axes will not be preserved.\n * Orphan bones - in direction from parent head\n * With only one non-orphan bone - to that one.",
            default = False,
            ) 
    import_mode = EnumProperty(
            name = "Import mode.",
            items = (('All','All','Import mesh and skeleton'),
                    ('Mesh','Mesh','Import only mesh'),
                    ('Skel','Skel','Import only skeleton'))
            )
    bDontInvertRoot= BoolProperty(
            name = "Don't invert root bone",
            description = " * Used by PSK and PSA.\n * Check it, if skeleton is badly oriented.",
            default = False,
            )
    bFilenameAsPrefix =  BoolProperty(
            name = "Prefix action names",
            description = "Use filename as prefix for action name.",
            default = False,
            )
    bActionsToTrack = BoolProperty(
            name = "All actions to NLA track",
            description = "Add all imported action to new NLAtrack. One by one.",
            default = False,
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
        if props.bDontInvertRoot:
            sub.label("", icon = 'ERROR')
            
        layout.prop(props, 'fBonesizeRatio')
        layout.prop(props, 'fBonesize')
        
    def draw_psa(self, context):
        props = context.scene.pskpsa_import
        layout = self.layout
        layout.prop(props,'bActionsToTrack')
        layout.prop(props,'bFilenameAsPrefix')
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

    filepath = StringProperty(
            subtype = 'FILE_PATH',
            )
    filter_glob = StringProperty(
            default = "*.psk;*.pskx",
            options = {'HIDDEN'},
            )
            
    def draw(self, context):
        self.draw_psk(context)
        # self.layout.prop(context.scene.pskpsa_import, 'bDontInvertRoot')
    # draw = ImportProps.draw_psk
    
    def execute(self, context):
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
        
        no_errors = pskimport( 
                        self.filepath,
                        context = context,
                        bImportmesh = bImportmesh, bImportbone = bImportbone,
                        fBonesize = props.fBonesize,
                        fBonesizeRatio = props.fBonesizeRatio,
                        bSpltiUVdata = props.bSpltiUVdata,
                        bReorientBones = props.bReorientBones,
                        bReorientDirectly = props.bReorientDirectly,
                        bDontInvertRoot = props.bDontInvertRoot,
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

    filepath = StringProperty(
            subtype = 'FILE_PATH',
            )
    filter_glob = StringProperty(
            default = "*.psa",
            options = {'HIDDEN'},
            )
            
    def draw(self, context):
        self.draw_psa(context)
        self.layout.prop(context.scene.pskpsa_import, 'bDontInvertRoot')
      
    def execute(self, context):
        props = context.scene.pskpsa_import
        psaimport(self.filepath,
            context = context,
            bFilenameAsPrefix = props.bFilenameAsPrefix, 
            bActionsToTrack = props.bActionsToTrack, 
            oArmature = blen_get_armature_from_selection(),
            bDontInvertRoot = props.bDontInvertRoot,
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

class PSKPSA_import_panel(bpy.types.Panel, ImportProps):
    bl_label = "PSK/PSA Import"
    bl_idname = "VIEW3D_PT_udk_import"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    
    # @classmethod
    # def poll(cls, context):
        # print(context.scene.get('pskpsa_import'),'poll')
        # context.scene.update_tag()
        # context.scene.update()
        # return context.scene.get('pskpsa_import') is not None
          
    def draw(self, context):
        props = context.scene.pskpsa_import
        if props is None:
            self.layout.label("??")
            return
        # return
        layout = self.layout
       
        # layout.label("Mesh and skeleton:")
        layout.operator(IMPORT_OT_psk.bl_idname, icon = 'MESH_DATA')
        self.draw_psk(context)
        # layout.prop(props, 'import_mode',expand = True)
        
        sub = layout.row()
        sub.operator(PSKPSA_OT_hide_unused_bones.bl_idname, icon = 'BONE_DATA')
        sub.enabled = (context.object is not
                       None) and (context.object.type == 'MESH' or context.object.type == 'ARMATURE')
                       
        layout.separator()
        layout.separator()
        # layout.label("Animation:", icon = 'ANIM')
        layout.operator(IMPORT_OT_psa.bl_idname, icon = 'ANIM')
        self.draw_psa(context)

    
def menu_func(self, context):
    self.layout.operator(IMPORT_OT_psk.bl_idname, text = "Skeleton Mesh (.psk)")
    self.layout.operator(IMPORT_OT_psa.bl_idname, text = "Skeleton Anim (.psa)")
    
    
def register():
    # print('register?')
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_import.append(menu_func)
    
    bpy.types.Scene.pskpsa_import = PointerProperty(type = PskImportOptions)
    
def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_import.remove(menu_func)
    
    del bpy.types.Scene.pskpsa_import
    
if __name__ == "__main__":
    register()

if __name__ == "io_import_scene_unreal_psa_psk_270_dev":
    import pskpsadev
