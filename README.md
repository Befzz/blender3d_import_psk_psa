Blender3D Import psk psa addon
========================
<ul>
<li>This is an heavily edited version of original blender plugin by Darknet / Optimus_P-Fat / Active_Trash / Sinsoft / flufy3d: https://en.blender.org/index.php/Extensions:2.6/Py/Scripts/Import-Export/Unreal_psk_psa (<a href="https://wiki.blender.org/index.php/Extensions:2.6/Py/Scripts/Import-Export/Unreal_psk_psa">old link</a>)
<li>Import mesh and skeleton from <b>.psk/.pskx</b></li>
<li>Import animation from <b>.psa</b></li>
<li>Game files can be exported to psk/psa by UModel: 
https://github.com/gildor2/UModel</li>
</ul>

<h5>Changes from original release</h5>
<ul>
<li>Blender 2.80+ support (experimental)</li>
<li>Fixed animation/skeleton import</li>
<li>Performance improvements</li>
<li>Panel UI updated</li>
<li>UI option: all actions to NLA track, one by one</li>
<li>UI option: mesh / skeleton or both import</li>
</ul>

<h3>Installation</h3>  

Disable original addon.

* Manually  

    Add .py file to the Blender's Add-ons search path:  
    * %APPDATA%\Blender Foundation\Blender\2.79\scripts\addons\  
    * %APPDATA%\Blender Foundation\Blender\2.80\scripts\addons\
  
* From Blender  
 
  2.79: File -> User preferences -> Add-ons -> Install Add-on from File...  
  2.80: Edit -> User preferences -> Add-ons -> Install Add-on from File...

<h3>Usage</h3>  

1. In 3DView, press **T** (Toggle Toolbar)  
2. Click **Misc.** tab.  

<table><tbody>
<tr><th> Panel in 3DView </th></tr>
<tr><td valign="top" align="center"><img src="https://github.com/Befzz/blender3d_import_psk_psa/blob/latest/imgs/panel.jpg"/></td>
</tr></tbody></table>
