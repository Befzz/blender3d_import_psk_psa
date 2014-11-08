blender3d_import_psk_psa
========================

import mesh, skeleton, animation from psk, psa files to blender3d

Tested on TERA(c) game files exported to psk/psa by UModel (https://github.com/gildor2/UModel)

<h4>Changes from original release</h4>

<ul>
<li>Animation import (semi-fixed, more test required)</li>
<li>Performance improvements</li>
<li>Code refractoring</li>
<li>Panel UI updated</li>
<li>UV textures names from material names</li>
<li>Option: prefix action name with filename</li>
<li>Option: mesh / bones or both import</li>
<li>Option: combined or separated UV maps</li>
</ul>
<h4>TODO</h4>
- [X] Support non ".psk/.psa" file extension. Checking by file header.
- [ ] Refract(nicefy,improve,style) psa import code
- [ ] Tests
- [ ] 100% correct of animation import

