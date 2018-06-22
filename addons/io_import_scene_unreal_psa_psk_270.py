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
    "name": "Import Unreal Skeleton Mesh (.psk)/Animation Set (.psa) (270) *latest",
    "author": "Darknet, flufy3d, camg188, befzz",
    "version": (2, 7, 1),
    "blender": (2, 64, 0),
    "location": "File > Import > Skeleton Mesh (.psk)/Animation Set (.psa) OR View3D > Tool Shelf (key T) > Misc. tab",
    "description": "Import Skeleleton Mesh / Animation Data",
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
- No Scale support. (no test material)
- Animation import updated (bone orientation now works)
"""

import bpy
import re
from mathutils import Vector, Matrix, Quaternion
from bpy.props import (FloatProperty,
                        StringProperty,
                        BoolProperty,
                        IntProperty,
                        EnumProperty,
                        PointerProperty )
                        
from struct import unpack, unpack_from
from bpy_extras.io_utils import unpack_list, unpack_face_list
import time

#DEV
# from mathutils import *
# from math import *

### Cross-version proxy functions 2.8 - 2.7 WiP
# bpy.app.version: (2,80,17)
v = bpy.app.version
if (v[0]*1000000 + v[1]*1000 + v[2]) >= 2080017:
  def scene_link_obj(obj):
    # return bpy.context.scene_collection.objects.link(obj)
    return bpy.context.collection.objects.link(obj)
    
  def util_obj_select(obj, action='SELECT'):
      if obj.name in bpy.context.view_layer.objects:
          return obj.select_set(action)
    
  def context_object_active(obj):
    bpy.context.view_layer.objects.active = obj
    
  def get_uv_layers(mesh_obj):
      return mesh_obj.uv_layers
else:
  def scene_link_obj(obj):
    return bpy.context.scene.objects.link(obj)
    
  def util_obj_select(obj, action='SELECT'):
    obj.select = (action == 'SELECT')
    
  def context_object_active(obj):
    bpy.context.scene.objects.active = obj

  def get_uv_layers(mesh_obj):
      return mesh_obj.uv_textures
del v

def utils_set_mode(mode):
    if bpy.ops.object.mode_set.poll():
        bpy.ops.object.mode_set(mode=mode, toggle = False)

# since names have type ANSICHAR(signed char) - using cp1251(or 'ASCII'?)
def util_bytes_to_str(in_bytes):
    return in_bytes.rstrip(b'\x00').decode(encoding='cp1252', errors='replace')

class class_psk_bone:
    bone_index = 0
    name = ""
    origmat = []
    scale = []
    parent = None
    parent_name = ""
    parent_index = 0
    __matrix_local_rot = None
    mat_world_rot = None
    __matrix_local = None
    mat_world = None
    trans = None
    rot = None
    gtrans = None
    grot = None
    
    orig_quat = None
    orig_loc = None
    
    children = None

        
def select_all(select):
    if select:
        actionString = 'SELECT'
    else:
        actionString = 'DESELECT'

    if bpy.ops.object.select_all.poll():
        bpy.ops.object.select_all(action=actionString)

    if bpy.ops.mesh.select_all.poll():
        bpy.ops.mesh.select_all(action=actionString)

    if bpy.ops.pose.select_all.poll():
        bpy.ops.pose.select_all(action=actionString)

def util_ui_show_msg(msg):
    bpy.ops.pskpsa.message('INVOKE_DEFAULT', message = msg)
        
        

#TODO check chunk flag?
def util_is_header_valid(filename, ftype, chunk_id, chunk_flag):
    '''Return True if chunk_id is a valid psk/psa (ftype) 'magick number'.'''
    if chunk_id != PSKPSA_FILE_HEADER[ftype]['chunk_id']:
        util_ui_show_msg(
            "The selected input file is not a " + ftype +
                        " file (header mismach)"
            "\nExpected: "+str(PSKPSA_FILE_HEADER[ftype]['chunk_id'])+
            "\nPresent: "+str(chunk_id)
        )    
        return False
    return True
    
    
def util_gen_name_part(filepath):
    '''strip path and extension from path'''
    return re.match(r'.*[/\\]([^/\\]+?)(\..{2,5})?$', filepath).group(1)
                
                
def direction_from_vec(vec_in, vec_out):
    vec = vec_in
    vec_ori = vec_out
    if abs(vec[0]) > abs(vec[1]):
        if abs(vec[0]) > abs(vec[2]):
            vec_ori.x = 1 if vec[0] >= 0 else -1
        else:
            vec_ori.z = 1 if vec[2] >= 0 else -1
    else:
        if abs(vec[1]) > abs(vec[2]):
            vec_ori.y = 1 if vec[1] >= 0 else -1
        else:
            vec_ori.z = 1 if vec[2] >= 0 else -1
            
PSKPSA_FILE_HEADER = {
    'psk':{'chunk_id':b'ACTRHEAD\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'},
    'psa':{'chunk_id':b'ANIMHEAD\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'}
}

class PskImporter:
    LEAF_AXIS_OWN_LOC = 0
    LEAF_FROM_PARENT = 1

  
  
    def calc_bone_rotation(self, psk_bone, blen_min, bonesize_auto, leaf_orientation = LEAF_AXIS_OWN_LOC):
        children = psk_bone.children

        vecy = Vector((0.0, 1.0, 0.0))
        quat = Quaternion((1.0, 0.0, 0.0, 0.0))
        vec_ori = Vector()
        
        # Leaf bone
        if len(children) == 0:
            # Single bone. ALONE.
            if psk_bone.parent == None:
                return (blen_min, quat)
                
            if leaf_orientation == self.LEAF_AXIS_OWN_LOC:
                direction_from_vec(psk_bone.trans.translation, vec_ori)
            elif leaf_orientation == self.LEAF_FROM_PARENT:
                vec_ori = psk_bone.oquat * psk_bone.trans.translation
              
            return (blen_min, vecy.rotation_difference(vec_ori))

        sumvec = Vector()
        for child in children:
            sumvec += (child.trans.translation)
        vec_to = sumvec / len(children)

        
        # print(vec, psk_bone.rot.to_euler(), psk_bone.name)
        # vec = vec_to * .7
        vec = vec_to
        direction_from_vec(vec, vec_ori)
        # print(psk_bone.name, vec_ori)
        
        # if psk_bone.name[:5] == 'thigh':
          # print(vec_to, vec, psk_bone.mat_world_rot.to_euler(), psk_bone.name)

        if bonesize_auto and vec.length > 0.01:
            blen = vec.length
        else:
            blen = blen_min

        return (blen, vecy.rotation_difference(vec_ori))

    
def pskimport(filepath, bImportmesh = True, bImportbone = True, bImportsingleuv = False, fBonesize = 2.0, 
        bDontInvertRoot = False,
        bReorientBones = False
        ):
    print(fBonesize,bReorientBones,bReorientBones)
    if not bImportbone and not bImportmesh:
        util_ui_show_msg("Nothing to do.\nSet something for import.")
        return False
    file_ext = 'psk'
    
    print ("-----------------------------------------------")
    print ("---------EXECUTING PSK PYTHON IMPORTER---------")
    print ("-----------------------------------------------")
    print (" Importing file:", filepath)

    #file may not exist
    try:
        pskfile = open(filepath,'rb')
    except IOError:
        util_ui_show_msg('Error while opening file for reading:\n  "'+filepath+'"')
        return False

    # logf.write('ChunkID: {0}\nTypeFlag: {1}\nDataSize: {2}\nDataCount: {3}\n'.format(
                # util_bytes_to_str(chunk_header_id),  chunk_header_type,
                # chunk_header_datasize, chunk_header_datacount))

    # using this instead of class to avoid "object.prop" lookup. 3x faster.
    chunk_header_id = None
    chunk_header_type = None
    chunk_header_datasize = None
    chunk_header_datacount = None
    #all binary Data of chunk for unpack (bytearray)
    chunk_data = None
    
    #=================================================
    #         VChunkHeader Struct
    # ChunkID|TypeFlag|DataSize|DataCount
    # 0      |1       |2       |3
    #=================================================
    # read("map"/"assign") a header and chunk data to local variables
    def read_chunk():
        nonlocal chunk_header_id,\
                 chunk_header_type,\
                 chunk_header_datasize,\
                 chunk_header_datacount,\
                 chunk_data
        #read header
        (chunk_header_id,
         chunk_header_type,
         chunk_header_datasize,
         chunk_header_datacount) = unpack('20s3i', pskfile.read(32))
        
        # print('HEADER',chunk_header_id, chunk_header_type, chunk_header_datasize, chunk_header_datacount)
        
        # read all chunk data
        chunk_data = pskfile.read(chunk_header_datacount * chunk_header_datasize)
    

    isPskx = False

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
        # printlog("New Mesh Data = " + mesh_data.name + "\n")
    #================================================================================================== 
    # General
    #================================================================================================== 
    read_chunk()
    
    # check file header
    if not util_is_header_valid(filepath, file_ext, chunk_header_id, chunk_header_type):
        return False

    #================================================================================================== 
    # Points (Vertices)
    #================================================================================================== 
    #read the PNTS0000 header ( VPoint )
    read_chunk()
    if bImportmesh:
        verts  = [None] * chunk_header_datacount
        
        for counter in range( chunk_header_datacount ):
            (vec_x, vec_y, vec_z) = unpack_from('3f', chunk_data, counter * chunk_header_datasize)
            verts[counter]  = (vec_x, vec_y, vec_z)
            
    #================================================================================================== 
    # Wedges (UV)
    #================================================================================================== 
    # https://github.com/gildor2/UModel/blob/master/Exporters/Psk.h
    # for struct of VVertex
    #
    #read the VTXW0000 header ( VVertex )
    read_chunk()
    
    if bImportmesh:
    
        if chunk_header_datacount > 65536:  # NumVerts
            print('Probably PSKX! %i Vertices.' % chunk_header_datacount)
            isPskx = True
            
        uv_material_indexes = []
        UVCoords = [None] * chunk_header_datacount
        #UVCoords record format = [pntIndx, U coord, v coord]
        # printlog("[pntIndx, U coord, v coord]\n");
        for counter in range( chunk_header_datacount ):
            (point_index,
             u, v,
             material_index) = unpack_from('=IffBxxx', chunk_data, counter * chunk_header_datasize )
             
            # print(point_index, u, v, material_index)
            UVCoords[counter] = (point_index, u, v)
            
            # printlog_line(point_index,u,v)
            if not material_index in uv_material_indexes:
                uv_material_indexes.append(material_index)
           
    #================================================================================================== 
    # Faces
    #================================================================================================== 
    #read the FACE0000 header
    read_chunk()
    if bImportmesh:
        #PSK FACE0000 fields: WdgIdx1|WdgIdx2|WdgIdx3|MatIdx|AuxMatIdx|SmthGrp
        #associate MatIdx to an image, associate SmthGrp to a material
        SGlist = []

        faces = [None] * chunk_header_datacount
        faceuv = [None] * chunk_header_datacount
        # facesmooth = []

        # smlist = []
        mat_groups = {}
        
        unpack_format = '=HHHBBI'
        if isPskx:
            unpack_format = '=IIIBBI'
            
        for counter in range(chunk_header_datacount):
            (pntIndxA, pntIndxB, pntIndxC,
             MatIndex, AuxMatIndex, SmoothingGroup
             ) = unpack_from(unpack_format, chunk_data, counter * chunk_header_datasize)
             
            # UVCoords is(point_index, u, v)
            #             0            1  2
            PNTSA = UVCoords[pntIndxC][0]
            PNTSB = UVCoords[pntIndxB][0]
            PNTSC = UVCoords[pntIndxA][0]
            #print(PNTSA, PNTSB, PNTSC) #face id vertex
 
            faces[counter] = (PNTSA, PNTSB, PNTSC, 0)

            uv = (
                ( UVCoords[pntIndxC][1], 1.0 - UVCoords[pntIndxC][2] ),
                ( UVCoords[pntIndxB][1], 1.0 - UVCoords[pntIndxB][2] ),
                ( UVCoords[pntIndxA][1], 1.0 - UVCoords[pntIndxA][2] )
            )
            
            faceuv[counter] = (uv, MatIndex, AuxMatIndex, SmoothingGroup)

            if not MatIndex in mat_groups:
                # print('mat:', MatIndex)
                mat_groups[MatIndex] = []
            mat_groups[MatIndex].append( uv )

            #collect a list of the smoothing groups
            # facesmooth.append(SmoothingGroup)

            # if not indata[5] in smlist:
                # print('SM:',indata[5])
                # smlist.append(indata[5])
                
            if SGlist.count(SmoothingGroup) == 0:
                SGlist.append(SmoothingGroup)
                # print("SmoothingGroups:", len(SGlist))
            #assign a material index to the face
            #Tmsh.faces[-1].materialIndex = SGlist.index(indata[5])
        # printlog("Using Materials to represent PSK Smoothing Groups...\n")
        
        # for mg in mat_groups:
            # print('mat_group,len:',mg,len(mat_groups[mg]))
    
    #================================================================================================== 
    # Materials
    #================================================================================================== 
    #read the MATT0000 header
    read_chunk()
    
    if bImportmesh:
        # print("-- Materials -- (index, name, faces)")
        materials = []
        
        for counter in range(chunk_header_datacount):

            (MaterialNameRaw,
             TextureIndex,
             PolyFlags,
             AuxMaterial,
             AuxFlags,
             LodBias,
             LodStyle ) = unpack_from('64s6i', chunk_data, chunk_header_datasize * counter)
            
            materialname = util_bytes_to_str( MaterialNameRaw )
            matdata = bpy.data.materials.get(materialname)
            
            if matdata is None:
                matdata = bpy.data.materials.new( materialname )
            # matdata = bpy.data.materials.new( materialname )
                
            materials.append( matdata)
            mesh_data.materials.append( matdata )
            # print(counter,materialname,TextureIndex)
            # if mat_groups.get(counter) is not None:
                # print("%i: %s" % (counter, materialname), len(mat_groups[counter]))

    #================================================================================================== 
    # Bones (VBone .. VJointPos )
    #================================================================================================== 
    #read the REFSKEL0 header
    #REFSKEL0 - Name|Flgs|NumChld|PrntIdx|Qw|Qx|Qy|Qz|LocX|LocY|LocZ|Lngth|XSize|YSize|ZSize
    read_chunk()

    psk_bones = []
    bni_dict = {}

    
    bone_by_parent_idx = [None] * chunk_header_datacount
    for counter in range( chunk_header_datacount ):
        
        (name_raw, flags, NumChildren, ParentIndex, #0 1 2 3
         quat_x, quat_y, quat_z, quat_w,            #4 5 6 7
         vec_x, vec_y, vec_z,                       #8 9 10
         joint_length,                              #11
         scale_x, scale_y, scale_z) = unpack_from('64s3i11f', chunk_data, chunk_header_datasize * counter)
    
        psk_bone = class_psk_bone()
        psk_bone.children = []
        
        bone_name = util_bytes_to_str(name_raw)

        psk_bone.name = bone_name
        psk_bone.bone_index = counter
        psk_bone.parent_index = ParentIndex
        
        # psk_bone.scale = (scale_x, scale_y, scale_z)

        bni_dict[psk_bone.name] = psk_bone.bone_index

        ### works
        mat_quat = Quaternion((quat_w, quat_x, quat_y, quat_z)).to_matrix().to_4x4()
        pos_vec = Vector((vec_x, vec_y, vec_z))
        ###

        psk_bone.trans = Matrix.Translation( pos_vec )
        # psk_bone.rot = mat_quat
        # mat_quat = Quaternion((quat_w, quat_x, quat_y, quat_z)).to_matrix().to_4x4()
        
        psk_bone.orig_quat = (quat_w, quat_x, quat_y, quat_z)
        psk_bone.orig_loc = (vec_x, vec_y, vec_z)

        psk_bone.oquat = Quaternion((quat_w, quat_x, quat_y, quat_z))

        if psk_bone.parent_index == 0 and psk_bone.bone_index == psk_bone.parent_index:
            if bDontInvertRoot:
                psk_bone.mat_world = mat_quat.copy()
                psk_bone.mat_world_rot = mat_quat.copy()
            else:
                psk_bone.mat_world = mat_quat.inverted()
                psk_bone.mat_world_rot = mat_quat.inverted()
                
            psk_bone.mat_world.translation = pos_vec

        psk_bones.append(psk_bone)
        #print(psk_bone.name, psk_bone.quater.to_quaternion())
        # print(counter, bone_name, ParentIndex)
        
    # root bone must have parent_index = 0 and selfindex = 0
    # TODO optimize math.
    for psk_bone in psk_bones:
            
        if psk_bone.parent_index == 0:
            if psk_bone.bone_index == 0:
                psk_bone.parent = None
                continue
        parent = psk_bones[psk_bone.parent_index]
        psk_bone.parent = parent
        psk_bone.parent_name = parent.name
        
        parent.children.append(psk_bone)
        
        # tr = Matrix.Translation(psk_bone.parent.mat_world.translation)
        # matrix_global = tr * ( psk_bone.parent.mat_world_rot * psk_bone.trans)

        psk_bone.mat_world =  ( parent.mat_world_rot * psk_bone.trans)
        psk_bone.mat_world.translation += parent.mat_world.translation
        matrix_global_rot = parent.mat_world_rot * psk_bone.oquat.conjugated().to_matrix().to_4x4()
        psk_bone.mat_world_rot = matrix_global_rot
        
   
    #psk_bones.sort( key=lambda bone: bone.parent_index)
    # print('-- Bones --')
    # print('Count: %i' % len(psk_bones))
    #================================================================================================
    # Blender armature
    #================================================================================================
    if fBonesize < 0.001:
        fBonesize = 0.001
        
    if bImportbone:
        armature_data = bpy.data.armatures.new(gen_names['armature_data'])
        armature_obj = bpy.data.objects.new(gen_names['armature_object'], armature_data)

        scene_link_obj(armature_obj)

        select_all(False)
        util_obj_select(armature_obj)
        context_object_active(armature_obj)
        
        utils_set_mode('EDIT')
        
        # TODO: options for axes and x_ray?
        armature_data.show_axes = False
        armature_data.draw_type = 'STICK'
        armature_obj.show_x_ray = True
        
        for psk_bone in psk_bones:
            edit_bone = armature_obj.data.edit_bones.new(psk_bone.name)

            armature_obj.data.edit_bones.active = edit_bone

            if psk_bone.parent is not None:
                edit_bone.parent = armature_obj.data.edit_bones[psk_bone.parent_name]
                
            
            if bReorientBones:
                (new_bone_size, quat_orient_diff) = PskImporter.calc_bone_rotation(PskImporter, psk_bone, fBonesize, True)
                post_quat = psk_bone.oquat.conjugated() * quat_orient_diff
            else:
                new_bone_size = fBonesize
                post_quat = psk_bone.oquat.conjugated()
            
            # only length is matter
            # edit_bone.tail = Vector((0.0, 0.0, new_bone_size))
            
            edit_bone.tail = Vector(( 0.0, new_bone_size, 0.0))
            # edit_bone.matrix = psk_bone.mat_world * quat_diff.to_matrix().to_4x4()
            edit_bone.matrix = psk_bone.mat_world * post_quat.to_matrix().to_4x4()
            
            
            # test_quat = Quaternion((Euler((-3,3,1.6))))
            
            
            #### FINAL
            # post_quat = psk_bone.oquat.conjugated() * quat_diff
            # edit_bone.matrix = psk_bone.mat_world * test_quat.to_matrix().to_4x4()
            # edit_bone["post_quat"] = test_quat
            #### 

            # edit_bone["post_quat"] = Quaternion((1,0,0,0))
            # edit_bone.matrix = psk_bone.mat_world* psk_bone.rot
     

            # if edit_bone.parent:
              # edit_bone.matrix = edit_bone.parent.matrix * psk_bone.trans * (psk_bone.oquat.conjugated().to_matrix().to_4x4())
              # edit_bone.matrix = edit_bone.parent.matrix * psk_bone.trans * (test_quat.to_matrix().to_4x4())
            # else:
              # edit_bone.matrix = psk_bone.oquat.to_matrix().to_4x4()
              
            # bindPose information for .psa
            edit_bone["post_quat"] = post_quat
            edit_bone["orig_quat"] = psk_bone.orig_quat
            edit_bone["orig_loc"] = psk_bone.orig_loc

    # utils_set_mode('OBJECT')

         
    #==================================================================================================
    #END BONE DATA BUILD
    #==================================================================================================
    '''
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
    '''
    #================================================================================================== 
    # Influences (Bone Weight)
    #================================================================================================== 
    #read the RAWW0000 header (VRawBoneInfluence)(Weight|PntIdx|BoneIdx)
    read_chunk()

    RWghts = [None] * chunk_header_datacount

    for counter in range(chunk_header_datacount):
        (Weight,
         PointIndex,
         BoneIndex ) = unpack_from('fii', chunk_data, chunk_header_datasize * counter)
         
        RWghts[counter] = (PointIndex, BoneIndex, Weight)
        
        #print("weight:", PointIndex, BoneIndex, Weight)

    # RWghts.sort( key=lambda wgh: wgh[0])
    
    # printlog("Vertex point and groups count = " + str(len(RWghts)) + "\n")
    # printlog("PntIdx\tBoneIdx\tWeight")
    
    # for vg in RWghts:
        # printlog(str(vg[0]) + "|" + str(vg[1]) + "|" + str(vg[2]) + "\n")

    """
    for x in range(len(Tmsh.faces)):
        for y in range(len(Tmsh.faces[x].v)):
            #find v in RWghts[n][0]
            findVal = Tmsh.faces[x].v[y].index
            n = 0
            while findVal != RWghts[n][0]:
                n = n + 1
            TmpCol = VtxCol[RWghts[n][1]]
            #check if a vertex has more than one influence
            if n != len(RWghts) - 1:
                if RWghts[n][0] == RWghts[n + 1][0]:
                    #if there is more than one influence, use the one with the greater influence
                    #for simplicity only 2 influences are checked, 2nd and 3rd influences are usually very small
                    if RWghts[n][2] < RWghts[n + 1][2]:
                        TmpCol = VtxCol[RWghts[n + 1][1]]
        Tmsh.faces[x].col.append(NMesh.Col(TmpCol[0], TmpCol[1], TmpCol[2], 0))
    """

    #================================================================================================== 
    # Building Mesh
    #================================================================================================== 
    if bImportmesh:
        mesh_data.vertices.add(len(verts))
        mesh_data.tessfaces.add(len(faces))
        mesh_data.vertices.foreach_set("co", unpack_list( verts ))
        mesh_data.tessfaces.foreach_set("vertices_raw", unpack_list( faces ))

        # for face in mesh_data.tessfaces:
            # .use_smooth is True or False - but facesmooth contains an int
            # TODO FIXME still incorrect
            # if facesmooth[face.index] > 0:
                # face.use_smooth = True

        utils_set_mode('OBJECT')

    #===================================================================================================
    # UV Setup
    #===================================================================================================
    if bImportmesh:
        if bImportsingleuv:
            get_uv_layers(mesh_data).new(name = "psk_uv_map_single")
            uvmap =  mesh_data.tessface_uv_textures[-1]
            # print("-- UV Single --\n" + uvmap.name)
            for face in mesh_data.tessfaces:
                face.material_index = faceuv[face.index][1]
                face_uv = faceuv[face.index][0]
                uvmap.data[face.index].uv1 = Vector((face_uv[0][0], face_uv[0][1]))
                uvmap.data[face.index].uv2 = Vector((face_uv[1][0], face_uv[1][1]))
                uvmap.data[face.index].uv3 = Vector((face_uv[2][0], face_uv[2][1]))
        else: #or make single UV map
            # print("-- UV Multi --")
            use_material_name = False
            if len(uv_material_indexes) == len(materials):
                use_material_name = True
            for i in range(len(uv_material_indexes)):
                
                if use_material_name:
                    uv_name = materials[i].name + ".uv"
                else:
                    uv_name = "psk_uv_multi_" + str(i)
                uv = get_uv_layers(mesh_data).new(name=uv_name)
                # print("%i: %s" % (i, uv.name))
                
            # creating different uv maps, if imported uv data have different uv_texture_id 
            _textcount = 0
            for uv in mesh_data.tessface_uv_textures:
                
                for face in mesh_data.tessfaces:
                    # faceuv is [] of (f_uv, MatIndex, AuxMatIndex, SmoothingGroup)
                    #                  0     1         2            3
                    # where f_uv is   ((u,v), (u,v), (u,v))
                    # if face index and texture index matches assign it
                    if faceuv[face.index][1] == _textcount:
                        mfaceuv = faceuv[face.index]
                        #assign material to face
                        face.material_index = faceuv[face.index][1]
                        
                        _uv1 = mfaceuv[0][0] #(0,0)
                        _uv2 = mfaceuv[0][1] #(0,0)
                        _uv3 = mfaceuv[0][2] #(0,0)
                        uv.data[face.index].uv1 = Vector((_uv1[0], _uv1[1])) #set them
                        uv.data[face.index].uv2 = Vector((_uv2[0], _uv2[1])) #set them
                        uv.data[face.index].uv3 = Vector((_uv3[0], _uv3[1])) #set them
                    else: #if not match zero them
                        uv.data[face.index].uv1 = Vector((0, 0)) #zero them 
                        uv.data[face.index].uv2 = Vector((0, 0)) #zero them 
                        uv.data[face.index].uv3 = Vector((0, 0)) #zero them 
                
                _textcount += 1
        #end if bImportsingleuv
        mesh_obj = bpy.data.objects.new(gen_names['mesh_object'], mesh_data)
    #===================================================================================================
    # Mesh Vertex Group bone weight
    #===================================================================================================
    if bImportmesh:
        #create bone vertex group #deal with bone id for index number
        for psk_bone in psk_bones:
            # group = mesh_obj.vertex_groups.new(bone.name)
            mesh_obj.vertex_groups.new(psk_bone.name)  
     
        for vgroup in mesh_obj.vertex_groups:
            # print(vgroup.name, ":", vgroup.index) 
            bone_index = bni_dict[vgroup.name]
            for vgp in RWghts:
                # vgp: 0, 1, 2 (vertexId, bone_index, weight)
                if vgp[1] == bone_index:
                    vgroup.add((vgp[0],), vgp[2], 'ADD')

        mesh_data.update()
        
        # bpy.context.scene.objects.link(mesh_obj) 
        scene_link_obj(mesh_obj)
        bpy.context.scene.update()

        select_all(False)
        #mesh_obj.select = True
        # bpy.context.scene.objects.active = mesh_obj
        context_object_active(mesh_obj)
    
    
        if bImportbone:
            bone_group_unused = armature_obj.pose.bone_groups.new(
                "Unused bones")
            bone_group_unused.color_set = 'THEME04'

            bone_group_nochild = armature_obj.pose.bone_groups.new(
                "No children")
            bone_group_nochild.color_set = 'THEME03'

            armature_data.show_group_colors = True

            for pbone in armature_obj.pose.bones:
                if mesh_obj.vertex_groups.find(pbone.name) == -1:
                    pbone.bone_group = bone_group_unused
                else:
                    if len(pbone.children) == 0:
                        pbone.bone_group = bone_group_nochild
                        
            # armature_obj.select = True
            # select_all(False)
            util_obj_select(armature_obj)
            # parenting mesh to armature object
            mesh_obj.parent = armature_obj
            mesh_obj.parent_type = 'OBJECT'
            # add armature modifier
            blender_modifier = mesh_obj.modifiers.new( armature_obj.data.name, type='ARMATURE')
            blender_modifier.show_expanded = False
            blender_modifier.use_vertex_groups = True
            blender_modifier.use_bone_envelopes = False
            blender_modifier.object = armature_obj
            
            # utils_set_mode('OBJECT')
            select_all(False)
            util_obj_select(armature_obj)
            context_object_active(armature_obj)
            # bpy.context.scene.update()
        
    utils_set_mode('OBJECT')
    return True


class class_psa_bone:
    name = ""
    # Transform = None
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

    # prev_quat = None
    # global_pos = None
    # global_rot = None

    
def blen_get_armature_from_selection():
  armature_obj = None
  for obj in bpy.context.selected_objects:
      if obj.type == 'ARMATURE':
          armature_obj = obj
          break
          
  if armature_obj is None:  
    for obj in bpy.context.selected_objects:
      if obj.type == 'MESH':
        for modifier in obj.modifiers:
          if modifier.type == 'ARMATURE':
            armature_obj = modifier.object
            break
    
  return armature_obj
  
    
def psaimport(filepath, bFilenameAsPrefix = False, bActionsToTrack = False, oArmature = None, 
              first_frames = 0, bDontInvertRoot = False):
    """Import animation data from 'filepath' using 'oArmature'
    
    Args:
        first_frames:
            Import only first_frames from each action
          
        bActionsToTrack:
            Put all imported action in one NLAtrack.
          
        oArmature:
            Skeleton used to calculate keyframes
    """
    print ("-----------------------------------------------")
    print ("---------EXECUTING PSA PYTHON IMPORTER---------")
    print ("-----------------------------------------------")
    print ("Importing file: ", filepath)
    
    file_ext = 'psa'
    try:
        psafile = open(filepath, 'rb')
    except IOError:
        util_ui_show_msg('Error while opening file for reading:\n  "'+filepath+'"')
        return False
        
    armature_obj = oArmature
    
    if armature_obj is None:  
        armature_obj = blen_get_armature_from_selection()
        if armature_obj is None:
            util_ui_show_msg("No armature selected.")
            return False


    chunk_header_id = None
    chunk_header_type = None
    chunk_header_datasize = None
    chunk_header_datacount = None
    chunk_data = None

    def read_chunk():
        nonlocal chunk_header_id, chunk_header_type,\
                 chunk_header_datasize, chunk_header_datacount,\
                 chunk_data

        (chunk_header_id, chunk_header_type,
         chunk_header_datasize, chunk_header_datacount) = unpack('20s3i', psafile.read(32))
        
        chunk_data = psafile.read(chunk_header_datacount * chunk_header_datasize)
    #==============================================================================================
    # General Header
    #==============================================================================================
    read_chunk()
    
    if not util_is_header_valid(filepath, file_ext, chunk_header_id, chunk_header_type):
        return False
    
    #==============================================================================================
    # Bones (FNamedBoneBinary)
    #==============================================================================================
    read_chunk()
    
    #Bones Data
    BoneIndex2NamePairMap = [None] * chunk_header_datacount
    BoneNotFoundList = []
    BonesWithoutAnimation = []
    
    # IndexToBone 

    # printlog("Name\tFlgs\tNumChld\tPrntIdx\tQx\tQy\tQz\tQw\tLocX\tLocY\tLocZ\tLength\tXSize\tYSize\tZSize\n")

    nobonematch = True
    
    # for case insensetive comparison
    # key = lowered name
    # value = orignal name
    skeleton_bones_lowered = {}
    
    for blender_bone_name in armature_obj.data.bones.keys():
      skeleton_bones_lowered[blender_bone_name.lower()] = blender_bone_name
      
    for counter in range(chunk_header_datacount):
        
        # tPrntIdx is -1 for parent; and 0 for other; no more useful data
        indata = unpack_from('64s3i11f', chunk_data, chunk_header_datasize * counter)

        bonename = util_bytes_to_str(indata[0])
        # bonename = util_bytes_to_str(indata[0]).upper()
        
        # print(counter, indata[2], indata[3], bonename)

        # if bonename in armature_obj.data.bones.keys():
        lowered = bonename.lower()
        if lowered in skeleton_bones_lowered:
            # BoneIndex2NamePairMap[counter] = bonename
            
            # use a skeleton bone name 
            BoneIndex2NamePairMap[counter] = skeleton_bones_lowered[lowered]
            #print('find bone', bonename)
            nobonematch = False
            
        else:
            # print("Can't find the bone:", bonename)
            BoneNotFoundList.append(counter)
            
    
    if nobonematch:
        util_ui_show_msg('No bone was match!\nSkip import!')
        return False
        
    for blender_bone_name in armature_obj.data.bones.keys():
        if BoneIndex2NamePairMap.count(blender_bone_name) == 0:
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
    Action_List = [None] * chunk_header_datacount
    
    for counter in range(chunk_header_datacount):
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
        ) = unpack_from('64s64s4i3f3i', chunk_data, chunk_header_datasize * counter)
                       
        action_name = util_bytes_to_str( action_name_raw )
        group_name = util_bytes_to_str( group_name_raw  )

        Raw_Key_Nums += Totalbones * NumRawFrames
        Action_List[counter] = ( action_name, group_name, Totalbones, NumRawFrames)

    #==============================================================================================
    # Raw keys (VQuatAnimKey)
    #==============================================================================================
    read_chunk()
    
    if(Raw_Key_Nums != chunk_header_datacount):
        util_ui_show_msg(
                'Raw_Key_Nums Inconsistent.'
                '\nData count found: '+chunk_header_datacount+
                '\nRaw_Key_Nums:' + Raw_Key_Nums
                )
        return False

    Raw_Key_List = [None] * chunk_header_datacount
    
    for counter in range(chunk_header_datacount):
        ( vec_x,  vec_y,  vec_z,
         quat_x, quat_y, quat_z, quat_w,
         time_until_next
        ) = unpack_from('3f4f1f', chunk_data, chunk_header_datasize * counter)
        
        pos = Vector((vec_x, vec_y, vec_z))
        quat = Quaternion((quat_w, quat_x, quat_y, quat_z))
        
        Raw_Key_List[counter] = (pos, quat, time_until_next)
        
    # Raw_Key_List = tuple(Raw_Key_List)
    
    utils_set_mode('OBJECT')

    #build tmp pose bone tree
    psa_bones = {}
    for pose_bone in armature_obj.pose.bones:
        psa_bone = class_psa_bone()
        psa_bone.name = pose_bone.name
        psa_bone.Transform = pose_bone.matrix
        if pose_bone.parent != None:
            psa_bone.parent = psa_bones[pose_bone.parent.name]
        else:
            psa_bone.parent = None
        psa_bones[pose_bone.name] = psa_bone
        psa_bone.orig_quat = Quaternion(pose_bone.bone['orig_quat']);
        psa_bone.post_quat = Quaternion(pose_bone.bone['post_quat']);
        psa_bone.orig_loc = Vector(pose_bone.bone['orig_loc']);

    # prepare
    # PsaBonesToProcess = tuple([psa_bones[x] for x in BoneIndex2NamePairMap if x is not None])
    # print(PsaBonesToProcess)
    
    # index of current frame in raw input data
    raw_key_index = 0
    
    ###dev
    def scene_update():
        bpy.context.scene.update()
    
    # unbind meshes, that uses this armature
    # because scene.update() calculating its positions
    # but we don't need it - its a big waste of time(CPU)
    
    armature_modifiers = []
    # for obj in bpy.data.objects:
        # if obj.type != 'MESH':
            # continue
        
        # for modifier in obj.modifiers:
            # if modifier.type != 'ARMATURE':
                # continue
            # if modifier.object == armature_obj:
                # armature_modifiers.append(modifier)
                # modifier.object = None
    
    armature_children = []
    
    #unbind children (same purpose)
    # for child in armature_obj.children:
        # armature_children.append((child, child.parent_type, child.parent_bone))
        # child.parent = None
    
    # bpy.context.scene.objects.active = armature_obj
    context_object_active(armature_obj)
    
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
    
    
    def _action_update_fcurves():
        for _, pbone in psa_bones.items():
            if pbone.fcurve_quat_w is None:
                continue
            pbone.fcurve_quat_w.update()
            pbone.fcurve_quat_x.update()
            pbone.fcurve_quat_y.update()
            pbone.fcurve_quat_z.update()
            pbone.fcurve_loc_x.update()
            pbone.fcurve_loc_y.update()
            pbone.fcurve_loc_z.update()
            
    counter = 0
    
    for raw_action in Action_List:
        ref_time = time.process_time()
        
        Name = raw_action[0]
        Group = raw_action[1]

        if Group != 'None':
            Name = "(%s) %s" % (Group,Name)
        if bFilenameAsPrefix:
            Name = "(%s) %s" % (gen_name_part, Name)
            
        Totalbones = raw_action[2]
        NumRawFrames = raw_action[3]
        
        action = bpy.data.actions.new(name = Name)

        
        # force print usefull information to console(due to possible long execution)
        counter += 1
        print("Action {0:>3d}/{1:<3d} frames: {2:>4d} {3}".format(
                counter, len(Action_List), NumRawFrames, Name)
              )
        
        #create all fcurves(for all bones) for an action
        for pose_bone in armature_obj.pose.bones:
            if pose_bone.name in BonesWithoutAnimation:
                continue
            psa_bone = psa_bones[pose_bone.name]
            
            data_path = pose_bone.path_from_id("rotation_quaternion")
            psa_bone.fcurve_quat_w = action.fcurves.new(data_path, index=0)
            psa_bone.fcurve_quat_x = action.fcurves.new(data_path, index=1)
            psa_bone.fcurve_quat_y = action.fcurves.new(data_path, index=2)
            psa_bone.fcurve_quat_z = action.fcurves.new(data_path, index=3)
        
            data_path = pose_bone.path_from_id("location")
            psa_bone.fcurve_loc_x = action.fcurves.new(data_path, index=0)
            psa_bone.fcurve_loc_y = action.fcurves.new(data_path, index=1)
            psa_bone.fcurve_loc_z = action.fcurves.new(data_path, index=2)
            
            # 1 Pre-add keyframes! \0/
            # 2 Set data: keyframe_points[].co[0..1]
            # 3 "Apply" data: fcurve.update() (before switching actions)
            psa_bone.fcurve_quat_w.keyframe_points.add(NumRawFrames)
            psa_bone.fcurve_quat_x.keyframe_points.add(NumRawFrames)
            psa_bone.fcurve_quat_y.keyframe_points.add(NumRawFrames)
            psa_bone.fcurve_quat_z.keyframe_points.add(NumRawFrames)

            psa_bone.fcurve_loc_x.keyframe_points.add(NumRawFrames) 
            psa_bone.fcurve_loc_y.keyframe_points.add(NumRawFrames) 
            psa_bone.fcurve_loc_z.keyframe_points.add(NumRawFrames) 
            
        pose_bones = armature_obj.pose.bones
        
        fcurve_interpolation = 'LINEAR'
        
        if first_frames > 0:
            maxframes = first_frames
        else:
            maxframes = 99999999
            
        # for i in range(NumRawFrames):
        for i in range(0,min(maxframes, NumRawFrames)):
            # raw_key_index+= Totalbones * 5 #55
            for j in range(Totalbones):
                if j in BoneNotFoundList:
                    raw_key_index += 1
                    continue
                
                bName = BoneIndex2NamePairMap[j]
                pbone = psa_bones[bName]
                pose_bone = pose_bones[bName]
                
                p_pos = Raw_Key_List[raw_key_index][0]
                p_quat = Raw_Key_List[raw_key_index][1]
                
                ##### Worked with no bone rotation
                # quat = p_quat.conjugated() * pbone.orig_quat
                # loc = p_pos - pbone.orig_loc
                #####
                if bDontInvertRoot:
                    quat = (p_quat * pbone.post_quat).conjugated() * (pbone.orig_quat * pbone.post_quat)
                else:
                    quat = (p_quat.conjugated() * pbone.post_quat).conjugated() * (pbone.orig_quat * pbone.post_quat)
                    
                loc = pbone.post_quat.conjugated() * p_pos -  pbone.post_quat.conjugated() * pbone.orig_loc

                if pbone.parent:
                    ##### Correct
                    # orig_prot = pose_bone.bone.parent.matrix_local.to_3x3().to_quaternion()
                    # orig_rot = pose_bone.bone.matrix_local.to_3x3().to_quaternion()
                    # orig_rot = (orig_prot.conjugated() * orig_rot)
                    ######

                    #### FINAL
                    quat = (p_quat * pbone.post_quat).conjugated() * (pbone.orig_quat * pbone.post_quat)
                    loc = pbone.post_quat.conjugated() * p_pos -  pbone.post_quat.conjugated() * pbone.orig_loc
                    ####
                    
                pose_bone.rotation_quaternion = quat
                pose_bone.location = loc
                    # pose_bone.rotation_quaternion = orig_rot.conjugated()
                    # pose_bone.location = p_pos - (pose_bone.bone.matrix_local.translation - pose_bone.bone.parent.matrix_local.translation)
                
                ##### Works + post_quat (without location works)
                # quat = (p_quat * pbone.post_quat).conjugated() * (pbone.orig_quat * pbone.post_quat)
                # loc = pbone.post_quat.conjugated() * (p_pos - pbone.orig_loc)

                
                pbone.fcurve_quat_w.keyframe_points[i].co = i, quat.w
                pbone.fcurve_quat_x.keyframe_points[i].co = i, quat.x
                pbone.fcurve_quat_y.keyframe_points[i].co = i, quat.y
                pbone.fcurve_quat_z.keyframe_points[i].co = i, quat.z
                
                # pbone.fcurve_quat_w.keyframe_points[i].interpolation = fcurve_interpolation
                # pbone.fcurve_quat_x.keyframe_points[i].interpolation = fcurve_interpolation
                # pbone.fcurve_quat_y.keyframe_points[i].interpolation = fcurve_interpolation
                # pbone.fcurve_quat_z.keyframe_points[i].interpolation = fcurve_interpolation
                
                pbone.fcurve_loc_x.keyframe_points[i].co = i, loc.x
                pbone.fcurve_loc_y.keyframe_points[i].co = i, loc.y
                pbone.fcurve_loc_z.keyframe_points[i].co = i, loc.z
                
                # pbone.fcurve_loc_x.keyframe_points[i].interpolation = fcurve_interpolation
                # pbone.fcurve_loc_y.keyframe_points[i].interpolation = fcurve_interpolation
                # pbone.fcurve_loc_z.keyframe_points[i].interpolation = fcurve_interpolation
                
                # Old path. Slower.
                # pbone.fcurve_quat_w.keyframe_points.insert(i,quat.w,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # pbone.fcurve_quat_x.keyframe_points.insert(i,quat.x,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # pbone.fcurve_quat_y.keyframe_points.insert(i,quat.y,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # pbone.fcurve_quat_z.keyframe_points.insert(i,quat.z,{'NEEDED','FAST'}).interpolation = fcurve_interpolation

                # pbone.fcurve_loc_x.keyframe_points.insert(i,loc.x,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # pbone.fcurve_loc_y.keyframe_points.insert(i,loc.y,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                # pbone.fcurve_loc_z.keyframe_points.insert(i,loc.z,{'NEEDED','FAST'}).interpolation = fcurve_interpolation
                raw_key_index += 1
            
            # on first frame
            # break
        raw_key_index += (NumRawFrames-min(maxframes,NumRawFrames))*Totalbones
        _action_update_fcurves()

            
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
            
        print("Done: %f sec.\n" % (time.process_time() - ref_time))
        #break on first animation set
        # break
    
    # set to rest position or set to first imported action
    if not bActionsToTrack:
        # pass
        if not bpy.context.scene.is_nla_tweakmode:
            armature_obj.animation_data.action = first_action
    # scene_update()
    
    # bind meshes again (setup modifier)
    for modifier in armature_modifiers:
        modifier.object = armature_obj
        
    # add children 
    for child in armature_children:
        (obj, p_type, p_bone) = child
        obj.parent = armature_obj
        obj.parent_type = p_type
        obj.parent_bone = p_bone
        
    select_all(False)
    util_obj_select(armature_obj)
    context_object_active(armature_obj)
    bpy.context.scene.frame_set(0)

 
class MessageOperator(bpy.types.Operator):
    bl_idname = "pskpsa.message"
    bl_label = "PSA/PSK"
    bl_options = {'REGISTER', 'INTERNAL'}

    message = StringProperty(default='Message')
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
        
        return context.window_manager.invoke_props_dialog(self, width=100 + 6*maxlen)
      
    def cancel(self, context):
        # print('cancel')
        self.execute(self)
        
    def draw(self, context):
        layout = self.layout
        sub = layout.column()
        sub.label(self.line0, icon='ERROR')

        for line in self.lines:
            sub.label(line)

# def getInputFilenamepsk(self, filename, bImportmesh, bImportbone, bImportsingleuv, fBonesize, bReorientBones):
def getInputFilenamepsk(self,filepath, **kwargs):
    return pskimport(        filepath, **kwargs)
    # return pskimport(         filename, bImportmesh, bImportbone, bImportsingleuv, fBonesize, bReorientBones):

# def getInputFilenamepsa(self, filename, _bFilenameAsPrefix, _bActionsToTrack):
  # return psaimport(         filename, bFilenameAsPrefix=_bFilenameAsPrefix, 
    # bActionsToTrack=_bActionsToTrack, 
    # oArmature = blen_get_armature_from_selection())
    
def getInputFilenamepsa(self, filepath, **kwargs):
  return psaimport(           filepath, **kwargs)
    
 
#properties for panels, and Operator.
class PskImportSharedOptions():
    # bl_options = {}
    # debug_log = BoolProperty(
            # name="Debug log",
            # description="Write debug information and raw data to <filename>.txt (for test purposes)",
            # default=False,
            # )
    fBonesize = FloatProperty(
            name="Bone length",
            description="Constant length for all bones. From head to tail distance",
            default=5.0, min=0.01, max=10, step=0.1, precision=2,
            )
    bImportsingleuv = BoolProperty(
            name="Single UV texture",
            description="If checked, MatIndex from vertex UV data will be ignored.",
            default=False,
            )
    bReorientBones = BoolProperty(
            name="Reorient bones",
            description="Bones will be axis-aligned to children.",
            default=False,
            )
    import_mode = EnumProperty(
            name="Import mode.",
            items=(('All','All','Import mesh and skeleton'),
                    ('Mesh','Mesh','Import only mesh'),
                    ('Skel','Skel','Import only skeleton'))
            )
    bDontInvertRoot= BoolProperty(
            name="Don't invert root bone",
            description=" * Used by PSK and PSA.\n * Check it, if skeleton is badly oriented.",
            default=False,
            )
    def draw_shared(self,opts):
        layout = self.layout
        layout.prop(opts, 'import_mode', expand=True)
        layout.prop(opts, 'bReorientBones')
        # layout.prop(opts, 'bDontInvertRoot')
        layout.prop(opts, 'bImportsingleuv')
        layout.prop(opts, 'fBonesize')
   
class PskImportOptions(bpy.types.PropertyGroup, PskImportSharedOptions):
    pass

def blen_hide_unused(armature_obj, mesh_obj):
    def is_bone_useless(pbone):
        is_useless = True
        if len(pbone.children) == 0:
            if mesh_obj.vertex_groups.get(pbone.name) != None:
                is_useless = False
        else:
            for pbone_child in pbone.children:
                is_useless = is_bone_useless(pbone_child)
                if not is_useless:
                    is_useless = False
                    # break
        # print(pbone.name, is_useless)
        if is_useless:
            pbone.hide = True
        return is_useless

    # print(armature_obj.data.bones[0].name)
    is_bone_useless(armature_obj.data.bones[0])


class ARMATURE_HIDE_UNUSED(bpy.types.Operator):
    """Hide useless bones(no weights and no children)\n* Select mesh with armature modifier.\n* ALT + H to reveal (in pose mode(CTRL + TAB))"""
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
                    print(
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

        
class IMPORT_OT_psk(bpy.types.Operator, PskImportSharedOptions):
    """Import skeleton and/or mesh."""
    bl_idname = "import_scene.psk"
    bl_label = "Import PSK"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_options = {'UNDO'}

    filepath = StringProperty(
            subtype='FILE_PATH',
            )
    filter_glob = StringProperty(
            default="*.psk;*.pskx",
            options={'HIDDEN'},
            )
            
    def draw(self, context):
        opts = bpy.context.scene.psk_import
        self.draw_shared(opts)

    def execute(self, context):
        opts = bpy.context.scene.psk_import
        if opts.import_mode == 'Mesh':
            bImportmesh = True
            bImportbone = False
        elif opts.import_mode == 'Skel':
            bImportmesh = False
            bImportbone = True
        else:
            bImportmesh = True
            bImportbone = True
        
        no_errors = getInputFilenamepsk(self, 
                        self.filepath,
                        bImportmesh = bImportmesh, bImportbone = bImportbone,
                        fBonesize = opts.fBonesize,
                        bImportsingleuv = opts.bImportsingleuv,
                        bReorientBones = opts.bReorientBones,
                        bDontInvertRoot = opts.bDontInvertRoot
                        )
        if not no_errors:
            return {'CANCELLED'}
        else:
            return {'FINISHED'}
    
    def invoke(self, context, event):
        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}

class IMPORT_OT_psa(bpy.types.Operator):
    '''Load a skeleton animation from .psa\n * Selected armature will be used.'''
    bl_idname = "import_scene.psa"
    bl_label = "Import PSA"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"

    filepath = StringProperty(
            subtype='FILE_PATH',
            )
    filter_glob = StringProperty(
            default="*.psa",
            options={'HIDDEN'},
            )
            
    def draw(self, context):
        layout = self.layout
        layout.prop(context.scene,'psa_bActionsToTrack')
        layout.prop(context.scene,'psa_bFilenameAsPrefix')
        layout.separator()
        layout.prop(context.scene.psk_import, 'bDontInvertRoot')
      
    def execute(self, context):
        getInputFilenamepsa(self, self.filepath,
            bFilenameAsPrefix = context.scene.psa_bFilenameAsPrefix, 
            bActionsToTrack = context.scene.psa_bActionsToTrack, 
            oArmature = blen_get_armature_from_selection(),
            bDontInvertRoot = context.scene.psk_import.bDontInvertRoot)
        return {'FINISHED'}
    
    def invoke(self, context, event):
        if blen_get_armature_from_selection() is None:
            util_ui_show_msg('Select an armature.')
            return {'FINISHED'}
        wm = context.window_manager
        wm.fileselect_add(self)
        return {'RUNNING_MODAL'}

bpy.types.Scene.psa_bFilenameAsPrefix =  BoolProperty(
            name="Prefix action names",
            description="Use filename as prefix for action name.",
            default=False,
            )
bpy.types.Scene.psa_bActionsToTrack = BoolProperty(
            name="All actions to NLA track",
            description="Add all imported action to new NLAtrack. One by one.",
            default=False,
            )

class Panel_UDKImport(bpy.types.Panel, PskImportSharedOptions):
    bl_label = "PSK/PSA Import"
    bl_idname = "VIEW3D_PT_udk_import"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    
    @classmethod
    def poll(cls, context):
        return context.scene.get('psk_import') is not None
          
    def draw(self, context):
        opts = context.scene.psk_import
        if opts is None:
            self.layout.label("??")
            return
        # return
        layout = self.layout
       
        # layout.label("Mesh and skeleton:")
        layout.operator(IMPORT_OT_psk.bl_idname, icon='MESH_DATA')
        self.draw_shared(opts)
        # layout.prop(opts, 'import_mode',expand=True)
        
        sub = layout.row()
        sub.operator(ARMATURE_HIDE_UNUSED.bl_idname, icon='BONE_DATA')
        sub.enabled = (context.object !=
                       None) and (context.object.type == 'MESH' or context.object.type == 'ARMATURE')
                       
        layout.separator()
        layout.label("Animation:")
        layout.operator(IMPORT_OT_psa.bl_idname, icon='ANIM')
        IMPORT_OT_psa.draw(self,context)
        # layout.prop(context.scene,'psa_bActionsToTrack')
        # layout.prop(context.scene,'psa_bFilenameAsPrefix')
    
def menu_func(self, context):
    self.layout.operator(IMPORT_OT_psk.bl_idname, text="Skeleton Mesh (.psk)")
    self.layout.operator(IMPORT_OT_psa.bl_idname, text="Skeleton Anim (.psa)")
    
    
def register():
    # print('register?')
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_import.append(menu_func)
    
    bpy.types.Scene.psk_import = PointerProperty(type=PskImportOptions)
    # dev_test()
    
def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_import.remove(menu_func)
    
    del bpy.types.Scene.psk_import
    
if __name__ == "__main__":
    register()

if __name__ == "io_import_scene_unreal_psa_psk_270_dev":
    import pskpsadev
