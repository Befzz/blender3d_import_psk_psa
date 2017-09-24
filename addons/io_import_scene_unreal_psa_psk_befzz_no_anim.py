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
    "name": "Import Unreal Skeleton Mesh (.psk) (no_anim)",
    # "name": "Import Unreal Skeleton Mesh (.psk)/Animation Set (.psa)",
    "author": "Darknet, flufy3d, camg188, befzz",
    "version": (2, 6, 4),
    "blender": (2, 64, 0),
    "location": "File > Import > Skeleton Mesh (.psk)/Animation Set (.psa)",
    "description": "Import Skeleleton Mesh/Animation Data",
    "warning": "",
    "wiki_url": "http://wiki.blender.org/index.php/Extensions:2.5/Py/"
                "Scripts/Import-Export/Unreal_psk_psa",
    "category": "Import-Export",
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
Version': '2.6.*' edited by befzz
Github: https://github.com/Befzz/blender3d_import_psk_psa
- No Scale support.
"""

import bpy
import math
import re
from mathutils import Vector, Matrix, Quaternion, Euler
from bpy.props import (FloatProperty,
                       StringProperty,
                       BoolProperty,
                       CollectionProperty,
                       IntProperty,
                       EnumProperty,
                       PointerProperty)
from bpy.types import UIList
from struct import unpack, unpack_from
from bpy_extras.io_utils import unpack_list, unpack_face_list


def utils_set_mode(mode):
    if bpy.ops.object.mode_set.poll():
        bpy.ops.object.mode_set(mode=mode, toggle=False)

# since names have type ANSICHAR(signed char) - using cp1251(or 'ASCII'?)


def util_bytes_to_str(in_bytes):
    return in_bytes.rstrip(b'\x00').decode(encoding='cp1252', errors='replace')


class class_md5_bone:
    bone_index = 0
    name = ""
    # bindpos = []
    #bindmat = []
    origmat = []
    scale = []
    parent = None
    parent_name = ""
    parent_index = 0
    matrix_local_rot = None
    matrix_global_rot = None
    matrix_local = None
    matrix_global = None
    # dev
    quat_local = None
    quat_global = None
    vec_tail_axis = Vector((0, 0, 1))

    have_weights = False

    # def __init__(self):
    # self.vec_tail_axis = Vector((0,0,1))

    def dump(self):
        print("bone index: ", self.bone_index)
        print("name: ", self.name)
        # print ("bind position: ", self.bindpos)
        print("parent: ", self.parent)
        print("parent index: ", self.parent_index)
        # print ("blenderbone: ", self.blenderbone)


def blen_calc_bone_orient(md5_bones, md5_bone, blen_min, bonesize_auto):
    children = []
    for bone in md5_bones:
        if bone.bone_index == 0:
            continue
        if bone.parent_index == md5_bone.bone_index:
            children.append(bone)

    # print("%s %i" % (md5_bone.name, len(children)))
    # print(md5_bone.name, len(children))
    if len(children) == 0:
        # print('leaf',md5_bone.vec_tail_axis)
        
        # Single bone
        if md5_bone.parent == None:                        
            return md5_bone.matrix_global * (md5_bone.vec_tail_axis * blen_min)
          
        if bonesize_auto:
            return md5_bone.matrix_global * (md5_bone.parent.vec_tail_axis.normalized() * blen_min)
        else:
            return md5_bone.matrix_global * md5_bone.parent.vec_tail_axis
        # dev

    sumvec = Vector()
    for child in children:
        sumvec += (child.matrix_global.translation)
    vec_to = sumvec / len(children) - md5_bone.matrix_global.translation

    vec = md5_bone.quat_global.inverted() * vec_to

    # mrot = md5_bone.matrix_global.to_3x3()

    if bonesize_auto:
        blen = vec.length
    else:
        blen = blen_min

    vec_ori = Vector()
    if abs(vec[0]) > abs(vec[1]):
        if abs(vec[0]) > abs(vec[2]):
            vec_ori.x = blen if vec[0] >= 0 else -blen
        else:
            vec_ori.z = blen if vec[2] >= 0 else -blen
    else:
        if abs(vec[1]) > abs(vec[2]):
            vec_ori.y = blen if vec[1] >= 0 else -blen
        else:
            vec_ori.z = blen if vec[2] >= 0 else -blen

    # vec_ori = vec_ori.normalized() * (vec.length/2)
    # if md5_bone.name == 'RightUpLeg':
    #     print('vec_ori', vec_ori)
    md5_bone.vec_tail_axis = vec_ori
    # print('XXX',md5_bone.vec_tail_axis)
    return md5_bone.matrix_global * vec_ori
    # return md5_bone.matrix_global * Vector((0,bonesize,0))


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
    bpy.ops.error.message_popup('INVOKE_DEFAULT', message=msg)


PSKPSA_FILE_HEADER = {
    'psk': {'chunk_id': b'ACTRHEAD\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'},
    'psa': {'chunk_id': b'ANIMHEAD\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00\x00'}
}
# TODO check chunk flag?


def util_is_header_valid(filename, ftype, chunk_id, chunk_flag):
    if chunk_id != PSKPSA_FILE_HEADER[ftype]['chunk_id']:
        util_ui_show_msg(
            "The selected input file is not a " + ftype +
            " file (header mismach)"
            "\nExpected: " + str(PSKPSA_FILE_HEADER[ftype]['chunk_id']) +
            "\nPresent: " + str(chunk_id)
        )
        return False
    return True


def util_gen_name_part(filepath):
    '''strip path and extension from path'''
    return re.match(r'.*[/\\]([^/\\]+?)(\..{2,5})?$', filepath).group(1)


def pskimport(filepath, bImportmesh, bImportbone, bDebugLogPSK, bImportsingleuv, bonesize, bonesize_auto):
    if not bImportbone and not bImportmesh:
        util_ui_show_msg("Nothing to do.\nSet something for import.")
        return False

    print("--------------------------------------------------")
    print("---------SCRIPT EXECUTING PYTHON IMPORTER---------")
    print("--------------------------------------------------")
    print(" Debug Enabled:", bDebugLogPSK)
    print(" Importing file:", filepath)

    # file may not exist
    try:
        pskfile = open(filepath, 'rb')
    except IOError:
        util_ui_show_msg(
            'Error while opening file for reading:\n  "' + filepath + '"')
        return False

    if bDebugLogPSK:
        logpath = filepath + ".txt"
        print("logpath:", logpath)
        logf = open(logpath, 'w')

    def printlog(strdata):
        if bDebugLogPSK:
            logf.write(strdata)

    def printlog_line(*args):
        if not bDebugLogPSK:
            return
        logf.write('\t'.join(list(map(str, args))) + "\n")

    def printlog_header():
        if not bDebugLogPSK:
            return
        logf.write('ChunkID: {0}\nTypeFlag: {1}\nDataSize: {2}\nDataCount: {3}\n'.format(
            util_bytes_to_str(chunk_header_id),  chunk_header_type,
            chunk_header_datasize, chunk_header_datacount))

    # using this instead of class to avoid "object.prop" lookup. 3x faster.
    chunk_header_id = None
    chunk_header_type = None
    chunk_header_datasize = None
    chunk_header_datacount = None
    # all binary Data of chunk for unpack (bytearray)
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
        # read header
        (chunk_header_id,
         chunk_header_type,
         chunk_header_datasize,
         chunk_header_datacount) = unpack('20s3i', pskfile.read(32))

        # print('HEADER',chunk_header_id, chunk_header_type, chunk_header_datasize, chunk_header_datacount)

        # read all chunk data
        chunk_data = pskfile.read(
            chunk_header_datacount * chunk_header_datasize)

        printlog_header()

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
        printlog("New Mesh Data = " + mesh_data.name + "\n")
    #==================================================================================================
    # General
    #==================================================================================================
    read_chunk()

    # check file header
    if not util_is_header_valid(filepath, 'psk', chunk_header_id, chunk_header_type):
        return False

    utils_set_mode('OBJECT')
    #==================================================================================================
    # Points (Vertices)
    #==================================================================================================
    # read the PNTS0000 header ( VPoint )
    read_chunk()
    if bImportmesh:
        verts = [None] * chunk_header_datacount

        for counter in range(chunk_header_datacount):
            (vec_x, vec_y, vec_z) = unpack_from(
                '3f', chunk_data, counter * chunk_header_datasize)
            # verts[counter]  = (vec_y, vec_x, vec_z)
            # FFIX
            verts[counter] = (vec_x, vec_y, vec_z)
            # EFFIX
            printlog_line(vec_x, vec_y, vec_z)

    #==================================================================================================
    # Wedges (UV)
    #==================================================================================================
    # https://github.com/gildor2/UModel/blob/master/Exporters/Psk.h
    # for struct of VVertex
    #
    # read the VTXW0000 header ( VVertex )
    read_chunk()

    if bImportmesh:

        if chunk_header_datacount > 65536:  # NumVerts
            print('Probably PSKX! %i Vertices.' % chunk_header_datacount)
            isPskx = True

        uv_material_indexes = []
        UVCoords = [None] * chunk_header_datacount
        # UVCoords record format = [pntIndx, U coord, v coord]
        printlog("[pntIndx, U coord, v coord, material_index]\n")
        for counter in range(chunk_header_datacount):
            (point_index,
             u, v,
             material_index) = unpack_from('=IffBxxx', chunk_data, counter * chunk_header_datasize)
            # =IffBxxx
            # = mean: don't align automatically

            # print(point_index, u, v, material_index)
            UVCoords[counter] = (point_index, u, v, material_index)

            printlog_line(point_index, u, v, material_index)
            if not material_index in uv_material_indexes:
                uv_material_indexes.append(material_index)

    #==================================================================================================
    # Faces
    #==================================================================================================
    # read the FACE0000 header or FACE3200 (pskx)
    read_chunk()
    if bImportmesh:
        # PSK FACE0000 fields: WdgIdx1|WdgIdx2|WdgIdx3|MatIdx|AuxMatIdx|SmthGrp
        # associate MatIdx to an image, associate SmthGrp to a material
        SGlist = []

        faces = [None] * chunk_header_datacount
        faceuv = [None] * chunk_header_datacount
        facesmooth = []
        printlog("nWdgIdx1\tWdgIdx2\tWdgIdx3\tMatIdx\tAuxMatIdx\tSmthGrp \n")

        # smlist = []
        mat_groups = {}

        unpack_format = '=HHHBBI'
        if isPskx:
            unpack_format = '=IIIBBI'

        for counter in range(chunk_header_datacount):
            (pntIndxA, pntIndxB, pntIndxC,
             MatIndex, AuxMatIndex, SmoothingGroup
             ) = unpack_from(unpack_format, chunk_data, counter * chunk_header_datasize)

            printlog_line(pntIndxA, pntIndxB, pntIndxC,
                          MatIndex, AuxMatIndex, SmoothingGroup)

            # UVCoords is(point_index, u, v)
            #             0            1  2
            PNTSA = UVCoords[pntIndxC][0]
            PNTSB = UVCoords[pntIndxB][0]
            PNTSC = UVCoords[pntIndxA][0]
            # print(PNTSA, PNTSB, PNTSC) #face id vertex

            faces[counter] = (PNTSA, PNTSB, PNTSC, 0)

            uv = (
                (UVCoords[pntIndxC][1], 1.0 - UVCoords[pntIndxC][2]),
                (UVCoords[pntIndxB][1], 1.0 - UVCoords[pntIndxB][2]),
                (UVCoords[pntIndxA][1], 1.0 - UVCoords[pntIndxA][2])
            )

            faceuv[counter] = (uv, MatIndex, AuxMatIndex, SmoothingGroup)

            if not MatIndex in mat_groups:
                # print('mat:', MatIndex)
                mat_groups[MatIndex] = []
            mat_groups[MatIndex].append(uv)

            # collect a list of the smoothing groups
            facesmooth.append(SmoothingGroup)

            # if not indata[5] in smlist:
            # print('SM:',indata[5])
            # smlist.append(indata[5])

            if SGlist.count(SmoothingGroup) == 0:
                SGlist.append(SmoothingGroup)
                # print("smooth:", SmoothingGroup)
            # assign a material index to the face
            #Tmsh.faces[-1].materialIndex = SGlist.index(indata[5])
        printlog("Using Materials to represent PSK Smoothing Groups...\n")

        # for mg in mat_groups:
        # print('mat_group,len:',mg,len(mat_groups[mg]))

    #==================================================================================================
    # Materials
    #==================================================================================================
    # read the MATT0000 header
    read_chunk()

    if bImportmesh:
        print("-- Materials -- (index, name, faces)")
        # printlog(" - Not importing any material data now. PSKs are texture wrapped! \n")
        materials = []

        for counter in range(chunk_header_datacount):

            (MaterialNameRaw,
             TextureIndex,
             PolyFlags,
             AuxMaterial,
             AuxFlags,
             LodBias,
             LodStyle) = unpack_from('64s6i', chunk_data, chunk_header_datasize * counter)

            materialname = util_bytes_to_str(MaterialNameRaw)
            matdata = bpy.data.materials.get(materialname)
            if matdata is None:
                matdata = bpy.data.materials.new(materialname)
            materials.append(matdata)
            mesh_data.materials.append(matdata)
            if mat_groups.get(counter) is not None:
                print("%i: %s" % (counter, materialname),
                      len(mat_groups[counter]))

    #==================================================================================================
    # Bones (VBone .. VJointPos )
    #==================================================================================================
    # read the REFSKEL0 header
    #REFSKEL0 - Name|Flgs|NumChld|PrntIdx|Qw|Qx|Qy|Qz|LocX|LocY|LocZ|Lngth|XSize|YSize|ZSize
    read_chunk()

    md5_bones = []
    bni_dict = {}

    printlog("Name\tFlgs\tNumChld\tPrntIdx\tQx\tQy\tQz\tQw\tLocX\tLocY\tLocZ\tLngth\tXSize\tYSize\tZSize\n")

    for counter in range(chunk_header_datacount):

        (name_raw, flags, NumChildren, ParentIndex,  # 0 1 2 3
         quat_x, quat_y, quat_z, quat_w,  # 4 5 6 7
         vec_x, vec_y, vec_z,  # 8 9 10
         joint_length,  # 11
         scale_x, scale_y, scale_z) = unpack_from('64s3I11f', chunk_data, chunk_header_datasize * counter)

        md5_bone = class_md5_bone()

        # print('vec_tail_axis', md5_bone.vec_tail_axis)
        bone_name = util_bytes_to_str(name_raw)

        printlog_line(util_bytes_to_str(name_raw), flags, NumChildren, ParentIndex, quat_x, quat_y, quat_z, quat_w,
                      vec_x, vec_y, vec_z, joint_length, scale_x, scale_y, scale_z)

        md5_bone.name = bone_name
        md5_bone.bone_index = counter
        md5_bone.parent_index = ParentIndex
        md5_bone.scale = (scale_x, scale_y, scale_z)

        # dev
        md5_bone.quat_local = Quaternion((quat_w, quat_x, quat_y, quat_z))
        md5_bone.quat_global = Quaternion((quat_w, quat_x, quat_y, quat_z))

        # dev
        # print(bone_name," %.2f, %.2f, %.2f" % tuple(math.degrees(a) for a in md5_bone.quat_local.to_euler()))

        bni_dict[md5_bone.name] = md5_bone.bone_index

        quat_mat = Quaternion(
            (quat_w, quat_x, quat_y, quat_z)).to_matrix().to_4x4()
        pos_vec = Vector((vec_x, vec_y, vec_z))

        # dev
        # if bone_name[:5] == 'Dummy':
        # print(bone_name, Quaternion((quat_w, quat_x, quat_y, quat_z)).to_euler())

        matrix_local_rot = quat_mat

        # don't back-rotate root bone.
        if md5_bone.parent_index == 0 and md5_bone.bone_index == md5_bone.parent_index:
            matrix_local = Matrix.Translation(pos_vec) * matrix_local_rot
            # matrix_local_rot = Euler((0,0,0)).to_matrix().to_4x4()
        else:
            matrix_local = Matrix.Translation(
                pos_vec) * matrix_local_rot.inverted()

        # matrix_local = matrix_local_rot * Matrix.Translation( pos_vec )

        md5_bone.matrix_local_rot = matrix_local_rot
        md5_bone.matrix_global_rot = matrix_local_rot
        md5_bone.matrix_local = matrix_local
        md5_bone.matrix_global = matrix_local

        md5_bones.append(md5_bone)
        #print(md5_bone.name, md5_bone.quater.to_quaternion())
        #print(counter, bone_name, ParentIndex)

    # root bone must have parent_index = 0 and selfindex = 0

    for md5_bone in md5_bones:

        if md5_bone.parent_index == 0:
            if md5_bone.bone_index == 0:
                md5_bone.parent = None
                continue

        md5_bone.parent = md5_bones[md5_bone.parent_index]
        md5_bone.parent_name = md5_bone.parent.name

        # global position + rotation
        matrix_global = md5_bone.parent.matrix_global * md5_bone.matrix_local
        md5_bone.matrix_global = matrix_global

        # matrix_global_rot = md5_bone.parent.matrix_global_rot * md5_bone.matrix_local_rot
        # matrix_global_rot =  md5_bone.matrix_local_rot * md5_bone.parent.matrix_global_rot
        matrix_global_rot = md5_bone.parent.matrix_global_rot * md5_bone.matrix_local_rot
        md5_bone.matrix_global_rot = matrix_global_rot

        md5_bone.quat_global = md5_bone.parent.quat_global.copy()
        # md5_bone.quat_global.rotate(md5_bone.quat_local)
        md5_bone.quat_global = md5_bone.parent.quat_global * \
            md5_bone.quat_local.conjugated()

        # if md5_bone.name[-3:] == 'HIP':
        # print(md5_bone.name)
        # print(md5_bone.quat_global.to_euler())
        # print(md5_bone.quat_local.to_euler())
        # print("%.2f, %.2f, %.2f" % tuple(math.degrees(a) for a in md5_bone.quat_global.to_euler()))
        # print("%.2f, %.2f, %.2f" % tuple(math.degrees(a) for a in md5_bone.quat_local.to_euler()))
    #md5_bones.sort( key=lambda bone: bone.parent_index)
    print('-- Bones --')
    print('Count: %i' % len(md5_bones))
    #================================================================================================
    # Blender armature
    #================================================================================================

    # obj = None
    # for obj in bpy.context.scene.objects:
    # if type(obj.data) is bpy.types.Armature:
    # armObj = obj
    # break

    # force create new armature if need
    if bImportbone:
        armature_data = bpy.data.armatures.new(gen_names['armature_data'])
        armature_obj = bpy.data.objects.new(
            gen_names['armature_object'], armature_data)

        bpy.context.scene.objects.link(armature_obj)
        # bpy.ops.object.mode_set(mode='OBJECT')

        select_all(False)
        armature_obj.select = True

        # set current armature to edit the bone
        bpy.context.scene.objects.active = armature_obj

        # TODO: options for axes and x_ray?
        armature_data.show_axes = True
        armature_data.draw_type = 'STICK'

        armature_obj.show_x_ray = True

        # Go to edit mode for the bones
        utils_set_mode('EDIT')

        for md5_bone in md5_bones:
            edit_bone = armature_obj.data.edit_bones.new(md5_bone.name)
            edit_bone.use_connect = False
            edit_bone.use_inherit_rotation = True
            edit_bone.use_inherit_scale = True
            edit_bone.use_local_location = True
            armature_obj.data.edit_bones.active = edit_bone

            ##########################################################
            joint_vector = md5_bone.matrix_global * Vector()
            edit_bone.head = joint_vector

            vector_tail_end_up = md5_bone.matrix_global_rot * Vector((0, 1, 0))
            # vector_tail_end_dir = md5_bone.matrix_global_rot * Vector((1,0,0))
            # vector_tail_end_dir = Vector((1,1,1)) * md5_bone.matrix_global_rot
            vector_tail_end_up.normalize()

            # dev
            # acopy = md5_bone.matrix_global #* md5_bone.quat_global.to_matrix().to_4x4()

            vecc = blen_calc_bone_orient(
                md5_bones, md5_bone, bonesize, bonesize_auto)

            zvec = Vector((0, 1, 0))

            edit_bone['quat_local'] = md5_bone.quat_local
            edit_bone['base_rotation'] = zvec.rotation_difference(
                md5_bone.vec_tail_axis)
            # edit_bone['base_rotation'] = md5_bone.quat_local
            # if md5_bone.

            # vecc = (md5_bone.matrix_global) * Vector((0,bonesize,0))
            # vecc.rotate(md5_bone.quat_global)
            # if md5_bone.parent_name == 'ROOT':
            # print(md5_bone.quat_global.to_euler())
            # print( md5_bone.matrix_global)
            # print( acopy)
            # print( vecc)
            edit_bone.tail = vecc

            edit_bone.align_roll(vector_tail_end_up)
            # edit_bone.align_roll(vecc)
            ###########################################################

            if md5_bone.parent is not None:
                edit_bone.parent = armature_obj.data.edit_bones[md5_bone.parent_name]
            else:
                continue

            '''
            # dev
            edit_bone = armature_obj.data.edit_bones.new(md5_bone.name + "__2")
            edit_bone.use_connect = False
            edit_bone.use_inherit_rotation = True
            edit_bone.use_inherit_scale = True
            edit_bone.use_local_location = True
            armature_obj.data.edit_bones.active = edit_bone
            
            obj = []
            vecfix = blen_calc_bone_orient(md5_bones, md5_bone, bonesize, bonesize_auto)
            
            # vecc = md5_bone.matrix_global * Vector((0,bonesize,0))
            vecc = md5_bone.matrix_global.translation + obj[0]
            joint_vector.y += 40
            vecc.y += 40
            edit_bone.head = joint_vector
            edit_bone.tail = vecc
            
            if md5_bone.parent is not None:
                edit_bone.parent = armature_obj.data.edit_bones[md5_bone.parent_name]
            '''

    # bpy.context.scene.update()
    #==================================================================================================
    # END BONE DATA BUILD
    #==================================================================================================
    '''
    if bImportmesh:
        VtxCol = []
        bones_count = len(md5_bones)
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
    # read the RAWW0000 header (VRawBoneInfluence)(Weight|PntIdx|BoneIdx)
    read_chunk()

    RWghts = [None] * chunk_header_datacount

    for counter in range(chunk_header_datacount):
        (Weight,
         PointIndex,
         BoneIndex) = unpack_from('fii', chunk_data, chunk_header_datasize * counter)

        RWghts[counter] = (PointIndex, BoneIndex, Weight)

        md5_bones[BoneIndex].have_weights = True

        #print("weight:", PointIndex, BoneIndex, Weight)

    RWghts.sort(key=lambda wgh: wgh[0])
    printlog("Vertex point and groups count = " + str(len(RWghts)) + "\n")
    printlog("PntIdx\tBoneIdx\tWeight")
    for vg in RWghts:
        printlog(str(vg[0]) + "|" + str(vg[1]) + "|" + str(vg[2]) + "\n")

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
    if (bDebugLogPSK):
        logf.close()
    #==================================================================================================
    # Building Mesh
    #==================================================================================================
    if bImportmesh:
        mesh_data.vertices.add(len(verts))
        mesh_data.tessfaces.add(len(faces))
        mesh_data.vertices.foreach_set("co", unpack_list(verts))
        mesh_data.tessfaces.foreach_set("vertices_raw", unpack_list(faces))

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
            mesh_data.uv_textures.new(name="psk_uv_map_single")
            uvmap = mesh_data.tessface_uv_textures[-1]
            print("-- UV Single --\n" + uvmap.name)
            for face in mesh_data.tessfaces:
                face.material_index = faceuv[face.index][1]
                face_uv = faceuv[face.index][0]
                uvmap.data[face.index].uv1 = Vector(
                    (face_uv[0][0], face_uv[0][1]))
                uvmap.data[face.index].uv2 = Vector(
                    (face_uv[1][0], face_uv[1][1]))
                uvmap.data[face.index].uv3 = Vector(
                    (face_uv[2][0], face_uv[2][1]))
        else:  # or make single UV map
            print("-- UV Multi --")
            use_material_name = False
            if len(uv_material_indexes) == len(materials):
                use_material_name = True
            for i in range(len(uv_material_indexes)):

                if use_material_name:
                    uv_name = materials[i].name + ".uv"
                else:
                    uv_name = "psk_uv_multi_" + str(i)
                uv = mesh_data.uv_textures.new(name=uv_name)
                print("%i: %s" % (i, uv.name))

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
                        # assign material to face
                        face.material_index = faceuv[face.index][1]

                        _uv1 = mfaceuv[0][0]  # (0,0)
                        _uv2 = mfaceuv[0][1]  # (0,0)
                        _uv3 = mfaceuv[0][2]  # (0,0)
                        uv.data[face.index].uv1 = Vector(
                            (_uv1[0], _uv1[1]))  # set them
                        uv.data[face.index].uv2 = Vector(
                            (_uv2[0], _uv2[1]))  # set them
                        uv.data[face.index].uv3 = Vector(
                            (_uv3[0], _uv3[1]))  # set them
                    else:  # if not match zero them
                        uv.data[face.index].uv1 = Vector((0, 0))  # zero them
                        uv.data[face.index].uv2 = Vector((0, 0))  # zero them
                        uv.data[face.index].uv3 = Vector((0, 0))  # zero them

                _textcount += 1
        # end if bImportsingleuv
        mesh_obj = bpy.data.objects.new(gen_names['mesh_object'], mesh_data)
    #===================================================================================================
    # Mesh Vertex Group bone weight
    #===================================================================================================
    if bImportmesh:
        # create bone vertex group #deal with bone id for index number
        for md5_bone in md5_bones:
            # group = mesh_obj.vertex_groups.new(bone.name)
            if md5_bone.have_weights:
                mesh_obj.vertex_groups.new(md5_bone.name)

        for vgroup in mesh_obj.vertex_groups:
            # print(vgroup.name, ":", vgroup.index)
            bone_index = bni_dict[vgroup.name]
            for vgp in RWghts:
                # vgp: 0, 1, 2 (vertexId, bone_index, weight)
                if vgp[1] == bone_index:
                    vgroup.add((vgp[0],), vgp[2], 'ADD')

        mesh_data.update()

        bpy.context.scene.objects.link(mesh_obj)
        bpy.context.scene.update()

        select_all(False)
        #mesh_obj.select = True
        bpy.context.scene.objects.active = mesh_obj

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

            armature_obj.select = True
            # parenting mesh to armature object
            mesh_obj.parent = armature_obj
            mesh_obj.parent_type = 'OBJECT'
            # add armature modifier
            blender_modifier = mesh_obj.modifiers.new(
                armature_obj.data.name, type='ARMATURE')
            blender_modifier.show_expanded = False
            blender_modifier.use_vertex_groups = True
            blender_modifier.use_bone_envelopes = False
            blender_modifier.object = armature_obj

    utils_set_mode('OBJECT')
    return True
#End of def pskimport#########################


class class_psa_bone:
    name = ""
    Transform = None
    parent = None
    fcurve_loc_x = None
    fcurve_loc_y = None
    fcurve_loc_z = None
    fcurve_quat_x = None
    fcurve_quat_y = None
    fcurve_quat_z = None
    fcurve_quat_w = None
    prev_quat = None
    # quat_local = None


def psaimport(filepath, context, bFilenameAsPrefix=False, bActionsToTrack=False, bArmatureSelected = False, armatureList = [], armatureListIdx = 0):
    print("--------------------------------------------------")
    print("---------SCRIPT EXECUTING PYTHON IMPORTER---------")
    print("--------------------------------------------------")
    print("Importing file: ", filepath)

    try:
        psafile = open(filepath, 'rb')
    except IOError:
        util_ui_show_msg(
            'Error while opening file for reading:\n  "' + filepath + '"')
        return False

    debug = True
    if (debug):
        logpath = filepath + ".txt"
        print("logpath:", logpath)
        logf = open(logpath, 'w')

    def printlog(strdata):
        if not debug:
            return
        logf.write(strdata)

    def printlogplus(name, data):
        if not debug:
            return

        logf.write(str(name) + '\n')
        if isinstance(data, bytes):
            # logf.write(str(bytes.decode(data).strip(bytes.decode(b'\x00'))))
            logf.write(util_bytes_to_str(data))
        else:
            logf.write(str(data))
        logf.write('\n')

    def write_log_plus(*args):
        if not debug:
            return
        for arg in args:
            if isinstance(arg, bytes):
                logf.write(util_bytes_to_str(arg) + '\t')
            else:
                logf.write(str(arg) + '\t')
        logf.write('\n')

    def write_log_plus_headers():
        write_log_plus(
            'ChunkID ',  chunk_header_id,
            'TypeFlag ', chunk_header_type,
            'DataSize ', chunk_header_datasize,
            'DataCount ', chunk_header_datacount)

    armature_obj = None

    if bArmatureSelected:
        # use selected armature
        if armatureList:
            armature_name = armatureList[armatureListIdx].name
            armature_obj = bpy.data.objects.get(armature_name)
            if armature_obj is None:
                util_ui_show_msg(
                    "Selected armature not found: " + armature_name)
                return False
    else:
        # use first armature
        for obj in bpy.data.objects:
            if obj.type == 'ARMATURE':
                armature_obj = obj
                break

    if armature_obj is None:
        util_ui_show_msg(
            "No armatures found.\nImport armature from psk file first.")
        if(debug):
            logf.close()
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

        chunk_data = psafile.read(
            chunk_header_datacount * chunk_header_datasize)
        write_log_plus_headers()
    #==============================================================================================
    # General Header
    #==============================================================================================
    read_chunk()

    if not util_is_header_valid(filepath, 'psa', chunk_header_id, chunk_header_type):
        if(debug):
            logf.close()
        return False

    #==============================================================================================
    # Bones (FNamedBoneBinary)
    #==============================================================================================
    read_chunk()

    # Bones Data
    BoneIndex2NamePairMap = [None] * chunk_header_datacount
    BoneNotFoundList = []
    BonesWithoutAnimation = []

    printlog("Name\tFlgs\tNumChld\tPrntIdx\tQx\tQy\tQz\tQw\tLocX\tLocY\tLocZ\tLength\tXSize\tYSize\tZSize\n")

    nobonematch = True

    for counter in range(chunk_header_datacount):
        indata = unpack_from('64s3i11f', chunk_data,
                             chunk_header_datasize * counter)

        bonename = util_bytes_to_str(indata[0])
        if bonename in armature_obj.data.bones.keys():
            BoneIndex2NamePairMap[counter] = bonename
            #print('find bone', bonename)
            nobonematch = False
        else:
            print('Can not find the bone:', bonename)
            BoneNotFoundList.append(counter)

    if nobonematch:
        util_ui_show_msg('No bone was match!\nSkip import!')
        if(debug):
            logf.close()
        return False

    for blender_bone_name in armature_obj.data.bones.keys():
        if BoneIndex2NamePairMap.count(blender_bone_name) == 0:
            BonesWithoutAnimation.append(blender_bone_name)
            print('Bone without animation frames:', blender_bone_name)

    #==============================================================================================
    # Animations (AniminfoBinary)
    #==============================================================================================
    read_chunk()

    Raw_Key_Nums = 0
    Action_List = [None] * chunk_header_datacount

    for counter in range(chunk_header_datacount):
        (action_name_raw,  # 0
         group_name_raw,  # 1
         Totalbones,  # 2
         RootInclude,  # 3
         KeyCompressionStyle,  # 4
         KeyQuotum,  # 5
         KeyReduction,  # 6
         TrackTime,  # 7
         AnimRate,  # 8
         StartBone,  # 9
         FirstRawFrame,  # 10
         NumRawFrames  # 11
         ) = unpack_from('64s64s4i3f3i', chunk_data, chunk_header_datasize * counter)

        write_log_plus('Name',        action_name_raw,
                       'Group',       group_name_raw,
                       'totalbones',  Totalbones,
                       'NumRawFrames', NumRawFrames
                       )

        action_name = util_bytes_to_str(action_name_raw)
        group_name = util_bytes_to_str(group_name_raw)

        Raw_Key_Nums += Totalbones * NumRawFrames
        Action_List[counter] = (action_name, group_name,
                                Totalbones, NumRawFrames)

    #==============================================================================================
    # Raw keys (VQuatAnimKey)
    #==============================================================================================
    read_chunk()

    if(Raw_Key_Nums != chunk_header_datacount):
        util_ui_show_msg(
            'Raw_Key_Nums Inconsistent.'
            '\nData count found: ' + chunk_header_datacount +
            '\nRaw_Key_Nums:' + Raw_Key_Nums
        )
        if(debug):
            logf.close()
        return False

    Raw_Key_List = [None] * chunk_header_datacount

    for counter in range(chunk_header_datacount):
        (vec_x,  vec_y,  vec_z,
         quat_x, quat_y, quat_z, quat_w,
         time_until_next
         ) = unpack_from('3f4f1f', chunk_data, chunk_header_datasize * counter)

        pos = Vector((vec_x, vec_y, vec_z))
        quat = Quaternion((quat_w, quat_x, quat_y, quat_z))

        Raw_Key_List[counter] = (pos, quat, time_until_next)

    utils_set_mode('OBJECT')

    # build tmp pose bone tree
    psa_bones = {}
    for bone in armature_obj.pose.bones:
        psa_bone = class_psa_bone()
        psa_bone.name = bone.name
        psa_bone.Transform = bone.matrix
        if bone.parent != None:
            psa_bone.parent = psa_bones[bone.parent.name]
        else:
            psa_bone.parent = None
        psa_bones[bone.name] = psa_bone

    print('Calculating animation:')

    # index of current frame in raw input data
    raw_key_index = 0

    # dev
    def scene_update():
        bpy.context.scene.update()

    # unbind meshes, that uses this armature
    # because scene.update() calculating its positions
    # but we don't need it - its a big waste of time(CPU)

    armature_modifiers = []
    for obj in bpy.data.objects:
        if obj.type != 'MESH':
            continue

        for modifier in obj.modifiers:
            if modifier.type != 'ARMATURE':
                continue
            if modifier.object == armature_obj:
                armature_modifiers.append(modifier)
                modifier.object = None

    armature_children = []
    # unbind children (same purpose)
    for child in armature_obj.children:
        armature_children.append((child, child.parent_type, child.parent_bone))
        child.parent = None

    # dev
    # for pose_bone in armature_obj.pose.bones:
        # pose_bone.bone.use_inherit_rotation = False
        # pose_bone.bone.use_inherit_scale = False

    bpy.context.scene.objects.active = armature_obj
    # armature_obj.hide = True
    # scene_update()

    ##########################################################
    mat_pose_rot_fix = Matrix.Rotation(-math.pi / 2,
                                       4, 'Z') * Matrix.Rotation(-math.pi / 2, 4, 'Y')
    ##########################################################

    counter = 0
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

    for raw_action in Action_List:
        counter += 1

        if counter != 4:
            continue

        Name = raw_action[0]
        Group = raw_action[1]

        if Group != 'None':
            Name = "(%s) %s" % (Group, Name)
        if bFilenameAsPrefix:
            Name = "(%s) %s" % (gen_name_part, Name)
        Totalbones = raw_action[2]
        NumRawFrames = raw_action[3]
        action = bpy.data.actions.new(name=Name)

        # force print usefull information to console(due to possible long execution)

        print("Action {0:>3d}/{1:<3d} frames: {2:>4d} {3}".format(
            counter, len(Action_List), NumRawFrames, Name)
        )

        # create all fcurves(for all bones) for frame
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

        pose_bones = armature_obj.pose.bones

        # for i in range(NumRawFrames):
        # first 5 frames(for testing)
        for i in range(0, 15):
            for j in range(Totalbones):
                if j not in BoneNotFoundList:
                    bName = BoneIndex2NamePairMap[j]
                    pbone = psa_bones[bName]
                    pose_bone = pose_bones[bName]

                    pos = Raw_Key_List[raw_key_index][0]
                    quat = Raw_Key_List[raw_key_index][1]
                    # print(pose_bone.name,pos)
                    # mat = Matrix()
                    quat_c = quat.conjugated()
                    # quat_c = quat
                    armature_obj.data.bones.active = pose_bone.bone

                    if pbone.parent != None:
                        base_rot = Quaternion(pose_bone.bone['base_rotation'])
                        quat_local = Quaternion(pose_bone.bone['quat_local'])

                        if pose_bone.name == 'RightArm':
                            print(pose_bone.bone.name,
                                  base_rot.to_euler(), quat.to_euler())
                        # matrix for calc's // calc from parent
                        # mat = Matrix.Translation(pos) * quat_c.to_matrix().to_4x4()
                        mat = (quat_c).to_matrix().to_4x4()
                        mat.translation = pos
                        mat = pbone.parent.Transform * mat

                        # matrix for posing
                        mat_view = pbone.parent.Transform * \
                            Matrix.Translation(pos)

                        # base_rot = armature_obj.data.bones[bNameq
                        # mat2 = (base_rot.inverted() * quat_c *
                        # base_rot).to_matrix().to_4x4()
                        # mat2.translation = pos
                        # mat2 = pbone.parent.Transform * mat2

                        # rot = quat_c * pbone.parent.Transform.to_quaternion()
                        # rot = rot.to_matrix().to_4x4()

                        # pose_bone.matrix = Matrix.Translation(mat_view.translation) \
                        # * rot
                        # pose_bone.matrix = mat2
                        origRot = quat_local
                        # pose_bone.location = pose_bone.bone.matrix_local.translation - pos

                        # eul = quat.to_euler()
                        # quatf = Euler((eul.x,eul.y,eul.z),'XYZ').to_quaternion()

                        # quatf = quat.rotation_difference(base_rot)
                        # quatf = origRot.rotation_difference(quat)
                        quatf = base_rot * quat * base_rot.inverted()

                        pose_bone.rotation_quaternion = quat
                        # pose_bone.rotation_quaternion = origRot.inverted() * quat * origRot
                        # pose_bone.rotation_quaternion = quat.inverted() * Quaternion()
                        # pose_bone.matrix = mat2 * base_rot.to_matrix().to_4x4()

                        # save mat for children calc's
                        pbone.Transform = mat
                    else:
                        # TODO fix needed?
                        mat = Matrix.Translation(
                            pos) * quat.to_matrix().to_4x4()
                        # pose_bone.matrix = mat * mat_pose_rot_fix
                        pbone.Transform = mat

                    loc = pose_bone.location
                    quat = pose_bone.rotation_quaternion

                    # update(calc) data (relative coordinates /location & rotation_quaternion/)
                    scene_update()

                    # possible fix for correct animation data(for interpolation)
                    if i != 0:
                        if pbone.prev_quat.dot(quat) < 0.0:
                            quat = -quat
                    pbone.prev_quat = quat.copy()

                    pbone.fcurve_quat_w.keyframe_points.insert(i, quat.w)
                    pbone.fcurve_quat_x.keyframe_points.insert(i, quat.x)
                    pbone.fcurve_quat_y.keyframe_points.insert(i, quat.y)
                    pbone.fcurve_quat_z.keyframe_points.insert(i, quat.z)

                    pbone.fcurve_loc_x.keyframe_points.insert(i, loc.x)
                    pbone.fcurve_loc_y.keyframe_points.insert(i, loc.y)
                    pbone.fcurve_loc_z.keyframe_points.insert(i, loc.z)

                raw_key_index += 1

        if bActionsToTrack:
            if nla_track_last_frame == 0:
                nla_stripes.new(Name, 0, action)
            else:
                nla_stripes.new(Name, nla_stripes[-1].frame_end, action)
            nla_track_last_frame += NumRawFrames
        elif is_first_action:
            first_action = action
            is_first_action = False

        # break on first animation set
        break

    # set to rest position or set to first imported action
    if not bActionsToTrack:
        # for pose_bone in armature_obj.pose.bones:
            # pose_bone.rotation_quaternion = (1,0,0,0)
            # pose_bone.location = (0,0,0)
        if not bpy.context.scene.is_nla_tweakmode:
            armature_obj.animation_data.action = first_action
    context.scene.frame_set(0)
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
    armature_obj.select = True
    bpy.context.scene.objects.active = armature_obj

    if(debug):
        logf.close()

    print('Done.')


class MessageOperator(bpy.types.Operator):
    bl_idname = "error.message_popup"
    bl_label = ""

    message = StringProperty(default='Message')
    lines = []

    def execute(self, context):
        self.lines = self.message.split("\n")
        maxlen = 0
        for line in self.lines:
            if len(line) > maxlen:
                maxlen = len(line)
            print(line)
        # self.lines.append("")
        wm = context.window_manager
        return wm.invoke_popup(self, width=30 + 6 * maxlen, height=400)

    def draw(self, context):
        layout = self.layout
        layout.label("[PSA/PSK Importer]", icon='ANIM')

        for line in self.lines:
            # row = self.layout.row(align=True)
            # row.alignment = 'LEFT'
            layout.label(line)


def getInputFilenamepsk(self, filename, bImportmesh, bImportbone, bDebugLogPSK, bImportsingleuv, bonesize, bonesize_auto):
    return pskimport(
        filename=filename,
        bImportmesh=bImportmesh,
        bImportbone=bImportbone,
        bDebugLogPSK=bDebugLogPSK,
        bImportsingleuv=bImportsingleuv,
        bonesize=bonesize,
        bonesize_auto=bonesize_auto)


def getInputFilenamepsa(self, filename, context, _bFilenameAsPrefix, _bActionsToTrack, bArmatureSelected, armatureList, armatureListIdx):
    return psaimport(filename,
        context=context,
        bFilenameAsPrefix=_bFilenameAsPrefix,
        bActionsToTrack=_bActionsToTrack,
        bArmatureSelected=bArmatureSelected,
        armatureList=armatureList,
        armatureListIdx=armatureListIdx)


class UDKImportArmaturePG(bpy.types.PropertyGroup):
    string = StringProperty()
    bones = StringProperty()
    have_animation = BoolProperty(default=False)

# properties for panels, and Operator.


class PskImportSharedOptions():
    bl_options = {}
    debug_log = BoolProperty(
        name="Debug log",
        description="Write debug information and raw data to <filename>.txt (for test purposes)",
        default=False,
    )
    bonesize = FloatProperty(
        name="Minimum Bone length",
        description="Constant length for all bones. From head to tail distance",
        default=3.0, min=0.1, max=10, step=0.1, precision=2,
    )
    single_uvtexture = BoolProperty(
        name="Single UV texture",
        description="If checked, MatIndex from vertex UV data will be ignored.",
        default=True,
    )
    import_mode = EnumProperty(
        name="Import mode.",
        items=(('All', 'All', 'Import mesh and skeleton'),
               ('Mesh', 'Mesh', 'Import only mesh'),
               ('Skel', 'Skel', 'Import only skeleton'))
    )


class PskImportOptions(bpy.types.PropertyGroup, PskImportSharedOptions):
    bonesize_auto = BoolProperty(
        name="Autosize bones",
        description="Bone size(length) determind by children[s]",
        default=True,
    )
    armature_list = CollectionProperty(type=UDKImportArmaturePG)
    armature_list_idx = IntProperty()
    armature_selected = BoolProperty(
        name="Armature Selected",
        description="Choose Armature to Import psa animation data",
        default=False)


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
    """Hide useless bones(no weights and no childrens)\n* Select mesh with armature modifier.\n* ALT + H to reveal (in pose mode(CTRL + TAB))"""
    bl_idname = "armature.hide_unused"
    bl_label = "Hide useless bones"
    bl_options = {'UNDO'}

    def execute(self, context):
        if context.object.type == 'MESH':
            for mod in context.object.modifiers:
                if mod.type == 'ARMATURE' and mod.object:
                    blen_hide_unused(
                        context.object.modifiers[0].object, context.object)
                else:

                    util_ui_show_msg(
                        "Can't find any situable Armature modifier for selected mesh.")
                    print(
                        "Can't find any situable Armature modifier for selected mesh.")
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
        layout = self.layout
        #layout.label("PSK Import", icon='ANIM')

        layout.prop(opts, 'import_mode', expand=True)
        sub = layout.row()
        sub.active = opts.import_mode != 'Skel'
        sub.prop(opts, 'single_uvtexture')
        sub = layout.row()
        sub.active = opts.import_mode != 'Mesh'
        sub.prop(opts, 'bonesize')
        sub = layout.row()
        sub.prop(opts, 'bonesize_auto')
        layout.prop(opts, 'debug_log')

    def execute(self, context):
        opts = bpy.context.scene.psk_import
        bImportbone = False
        bImportmesh = False
        if opts.import_mode == 'Mesh':
            bImportmesh = True
        elif opts.import_mode == 'Skel':
            bImportbone = True
        else:
            bImportbone = True
            bImportmesh = True

        no_errors = getInputFilenamepsk(self,
                                        self.filepath,
                                        bImportmesh, bImportbone, opts.debug_log,
                                        opts.single_uvtexture,
                                        bpy.context.scene.psk_import.bonesize,
                                        bpy.context.scene.psk_import.bonesize_auto
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
    '''Load a skeleton anim psa File'''
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
    bFilenameAsPrefix = BoolProperty(
        name="Prefix action names",
        description="Use filename as prefix for action names.",
        default=False,
    )
    bActionsToTrack = BoolProperty(
        name="All actions to NLA track",
        description="Add all imported actions to new NLAtrack. One by one.",
        default=False,
    )

    def execute(self, context):
        # getInputFilenamepsa(self, self.filepath, context,
                            # self.bFilenameAsPrefix, self.bActionsToTrack)
        return {'FINISHED'}

    def invoke(self, context, event):
        # wm = context.window_manager
        # wm.fileselect_add(self)
        return {'RUNNING_MODAL'}


class Panel_UDKImport(bpy.types.Panel):
    bl_label = "PSK/PSA Import"
    bl_idname = "OBJECT_PT_udk_import"
    bl_space_type = "VIEW_3D"
    bl_region_type = "TOOLS"
    #bl_options = {'HIDE_HEADER'}
    # filepath = StringProperty(
    # subtype='FILE_PATH',
    # )

    # @classmethod
    # def poll(cls, context):
    # return context.scene.get('psk_import') is not None

    def draw(self, context):
        opts = bpy.context.scene.psk_import
        layout = self.layout

        layout.label("Mesh and skeleton:")
        layout.operator(IMPORT_OT_psk.bl_idname, icon='MESH_DATA')
        layout.prop(opts, 'import_mode', expand=True)

        sub = layout.row()
        sub.operator(ARMATURE_HIDE_UNUSED.bl_idname, icon='BONE_DATA')
        sub.enabled = (context.object !=
                       None) and context.object.type == 'MESH'

        # layout.separator()
        # layout.label("Animation:")
        # layout.operator(IMPORT_OT_psa.bl_idname, icon='ANIM')

        # if opts.armature_selected:
            # split = layout.split(.75)
            # split.prop(opts, "armature_selected")
            # split.operator(OBJECT_OT_UDKImportArmature.bl_idname,
                           # text="", icon='FILE_REFRESH')
            # layout.template_list("OBJECT_UL_armatures", "", opts, "armature_list",
                                 # opts, "armature_list_idx", rows=5)
        # else:
            # layout.prop(opts, "armature_selected")


class OBJECT_UL_armatures(UIList):
    def draw_item(self, context, layout, data, item, icon, active_data, active_propname, index):
        split = layout.split(0.75)
        if item.have_animation:
            split.label(str(item.name), icon='ANIM')
        else:
            split.label("     " + str(item.name))
        #split.prop(item, "bones", text="", emboss=False, translate=False, icon='BONE_DATA')
        split.label(str(item.bones), icon='BONE_DATA')


class OBJECT_OT_PSAPath(bpy.types.Operator):
    """Select .psa file path to import for animation data"""
    bl_idname = "object.psapath"
    bl_label = "Import PSA Path"

    filepath = StringProperty(
        name="PSA File Path",
        description="Filepath used for importing the PSA file",
        maxlen=1024,
        default=""
    )
    filter_glob = StringProperty(
        default="*.psa",
        options={'HIDDEN'},
    )
    bFilenameAsPrefix = BoolProperty(
        name="Prefix action names",
        description="Use filename as prefix for action names.",
        default=False
    )

    def execute(self, context):
        getInputFilenamepsa(self, self.filepath, context)
        return {'FINISHED'}

    def invoke(self, context, event):
        bpy.context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}


class OBJECT_OT_UDKImportArmature(bpy.types.Operator):
    """Update armature list"""
    bl_idname = "object.udkimportarmature"
    bl_label = ""  # not visible in "<space>"

    def execute(self, context):
        my_objlist = bpy.context.scene.psk_import.armature_list
        objectl = []

        # clear list
        for _obj in my_objlist:
            my_objlist.remove(0)

        # fill list
        for obj in bpy.context.scene.objects:
            if obj.type == 'ARMATURE':
                list_item = my_objlist.add()
                list_item.name = obj.name
                list_item.bones = str(len(obj.data.bones))

                # if obj have assigned action or any NLA tracks, mark it
                if obj.animation_data\
                        and (obj.animation_data.action
                             or obj.animation_data.nla_tracks):
                    list_item.have_animation = True
                else:
                    list_item.have_animation = False

        return{'FINISHED'}


def menu_func(self, context):
    self.layout.operator(IMPORT_OT_psk.bl_idname, text="Skeleton Mesh (.psk)")
    # self.layout.operator(IMPORT_OT_psa.bl_idname, text="Skeleton Anim (.psa)")


def register():
    bpy.utils.register_module(__name__)
    bpy.types.INFO_MT_file_import.append(menu_func)
    bpy.types.Scene.psk_import = PointerProperty(
        type=PskImportOptions, name="Psk/Psa import options")


def unregister():
    bpy.utils.unregister_module(__name__)
    bpy.types.INFO_MT_file_import.remove(menu_func)
    del bpy.types.Scene.psk_import


if __name__ == "__main__":
    register()
