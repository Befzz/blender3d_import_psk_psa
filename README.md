Blender3D Import psk psa addon
========================
<ul>
<li>This is an heavily edited version of original blender plugin: https://wiki.blender.org/index.php/Extensions:2.6/Py/Scripts/Import-Export/Unreal_psk_psa
<li>Import mesh and skeleton from <b>.psk/.pskx</b></li>
<li>Import animation from <b>.psa</b></li>
<li>Game files can be exported to psk/psa by UModel: 
https://github.com/gildor2/UModel</li>
</ul>
<h4>Version 280</h4>

 * Blender 2.80+ support (experimental).
 * Animation + bone orientation fixed.
 * A part of code has been revised.


<h3>Versions</h3>

. | befzz | no_anim | 280
------------ | ------ | ------- | --
Animation | Bad | - | **Yes**
Bone orientation | Bad | Axis based | **Axis + Direct**
Bone autosize | - | Yes | **Yes**

<h5>Changes from original release</h5>
<ul>
<li>Fixed animation/skeleton import(280 version)</li>
<li>Performance improvements</li>
<li>Panel UI updated</li>
<li>UI option: prefix action name with filename</li>
<li>UI option: all actions to NLA track, one by one</li>
<li>UI option: mesh / bones or both import</li>
<li>UI option: combined or separated UV maps</li>
</ul>

<h3>Installation</h3>
Disable original addon.

* Manually  

    Add .py files to  
   *%APPDATA%\Blender Foundation\Blender\\* &nbsp;&nbsp; *VERSION* &nbsp;&nbsp; *\scripts\addons* 
  
* From blender  
 
  a) *Import from File* (File -> User preferences -> Add-ons -> Install from File...)  
  
  b) Or press space in 3D-view and enter "ad ." (ad, space, dot)

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
<td valign="top" align="center"><img alt="[Panel in 3D view]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/latest/imgs/panel.jpg"/></td>
<td valign="top" align="center"><img alt="[Panel in file selector]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/latest/imgs/psk_file_options.jpg"/></td>
</tr>
<tr><th colspan="2">Prefixed actions imported to one NLA track</th></tr>
<tr><td colspan="2" valign="top" align="center"><img alt="[Panel in 3D view]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/latest/imgs/nla_track.jpg"/>
</td></tr>
<tr><th colspan="2">Assembled models (pose from animation)</th></tr>
<tr><td align="center">TERA: The Exiled Realm of Arborea</td><td align="center">Alice. Madness Returns</td></tr>
<tr><td valign="top"><img alt="[TERA: The Exiled Realm of Arborea]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/latest/imgs/tera_test.jpg"/></td>
<td valign="top"><img alt="[Alice. Madness Returns]" src="https://github.com/Befzz/blender3d_import_psk_psa/blob/latest/imgs/alice_test.jpg"/></td>
</tr></tbody></table>
