Blender3D Import psk psa addon
========================

<h5>Versions</h5>

. | befzz | befzz no_anim
------------ | ------ | -------
Animation | **Yes** | No
Bone orientation | Bad | **Good**

<h5>Description</h5>
<ul>
<li>Mesh and skeleton from <b>.psk/.pskx</b></li>
<li>Animation from <b>.psa</b></li>
</ul>
Game files can be exported to psk/psa by UModel: 
https://github.com/gildor2/UModel
<h5>Changes from original release</h5>
<ul>
<li>Animation import (fixed partially)</li>
<li>Performance improvements</li>
<li>Panel UI updated</li>
<li>UV textures names from material names</li>
<li>Option: prefix action name with filename</li>
<li>Option: all actions to NLA track, one by one</li>
<li>Option: mesh / bones or both import</li>
<li>Option: combined or separated UV maps</li>
</ul>
<h5>Not supported</h5>
<ul>
<li>Bone scale</li>
</ul>
</ul>
<h5>Tested on(animation and mesh)</h5>
<ul>
<li>TERA: The Exiled Realm of Arborea</li>
<li>Alice. Madness Returns</li>
<li>Life is Strange - Episode 1</li>
<li>PUBG (Aug 2017) (.pskx mesh tested)</li>
</ul>

<table><tbody>
<tr><th> Panel in 3D view </th><th> Panel in file selector for psk </th></tr>
<tr>
<td valign="top" align="center"><img alt="[Panel in 3D view]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/master/imgs/panel.jpg"/></td>
<td valign="top" align="center"><img alt="[Panel in file selector]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/master/imgs/psk_file_options.jpg"/></td>
</tr>
<tr><th colspan="2">Prefixed actions imported in one NLA track</th></tr>
<tr><td colspan="2" valign="top" align="center"><img alt="[Panel in 3D view]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/master/imgs/nla_track.jpg"/>
</td></tr>
<tr><th colspan="2">Assembled models in pose from animation</th></tr>
<tr><td align="center">TERA: The Exiled Realm of Arborea</td><td align="center">Alice. Madness Returns</td></tr>
<tr><td valign="top"><img alt="[TERA: The Exiled Realm of Arborea]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/master/imgs/tera_test.jpg"/></td>
<td valign="top"><img alt="[Alice. Madness Returns]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/master/imgs/alice_test.jpg"/></td>
</tr></tbody></table>
