# import_trajectory.py
# ---------------------------------------------------------------------------
# Import a coupled Mem3DG + Cytosim trajectory into Blender as an animation.
#
#   * membrane : a mesh, one OBJ per frame   (frames/membrane_%04d.obj)
#   * fibers   : curves, one text file per frame (fibers/fibers_%04d.txt)
#   * points   : singles/couples, one text file per frame (points_%04d.txt),
#                drawn as vertex-instanced spheres coloured by binding state
#
# All three sequences share a frame index: membrane_000K.obj, fibers_000K.txt
# and points_000K.txt are
# the SAME simulated instant.  Coordinates are raw simulation micrometres in a
# shared world frame -- see COORDINATES.md.  One global WORLD_SCALE is applied
# so the whole scene is a convenient size in Blender; it does NOT change the
# relative geometry.
#
# USAGE
#   1. edit the CONFIG block below (paths, scale)
#   2. Blender:  Scripting workspace -> Open -> this file -> Run
#      or headless:  blender --background --python import_trajectory.py
#
# The animation is driven by a frame-change handler: as Blender's timeline
# advances, the membrane mesh is swapped and the fiber curves and point clouds
# are rebuilt for that frame.  Scrub the timeline or press Play.
# ---------------------------------------------------------------------------

import bpy
import os
import glob
import bmesh

# ============================ CONFIG =======================================
# Directories holding the two sequences (absolute paths recommended).
#
# When the frames come from coupled_driver both sequences land in the SAME
# directory, so these two are equal -- that is the default below, pointing at
# the bundled test dataset.  They are kept separate so you can also mix a
# membrane run with fibers exported from a standalone Cytosim trajectory via
# the `export_fibers` tool.
# surface nucleation dataset (membrane + fibers + singles/couples).  Swap to
# ".../test-membrane-response/frames-remesh" for the remeshed protrusion run.
# _TEST = "/Users/hakenuwa/Desktop/Mem3DG/tests/test-membrane-response/frames-remesh"
# _TEST = "/Users/hakenuwa/Desktop/Mem3DG/tests/test-multi-membrane/frames-approach"
# _TEST = "/Users/hakenuwa/Desktop/Mem3DG/tests/test-multi-membrane/frames-inflate"
_TEST = "/Users/hakenuwa/Desktop/Mem3DG/tests/test-nucleation-anchor/frames-steric"

MEMBRANE_DIR = _TEST     # membrane_%04d.obj
FIBER_DIR    = _TEST     # fibers_%04d.txt

# One scale for the whole scene.  Simulation units are ~micrometres; values are
# small (order 1), so 1.0 is usually fine.  This is the ONLY transform applied;
# keep it identical for membrane and fibers so they stay registered.
WORLD_SCALE  = 1.0

POINT_DIR    = _TEST     # points_%04d.txt  (singles / couples)

# Fiber curve appearance.
FIBER_BEVEL_FROM_RADIUS = True   # use each fiber's radius as bevel depth
FIBER_BEVEL_FALLBACK    = 0.01   # bevel depth if radius is 0 / disabled

# Singles/couples appearance.  Cytosim gives these no physical radius, so the
# on-screen size is purely a display choice.  Per-class overrides are keyed by
# the class name from the file ("single surfaceAnchor free", or just the
# property name "surfaceAnchor" to catch every state of that class).
POINT_RADIUS          = 0.012
POINT_RADIUS_OVERRIDE = {}
# Colour per binding state, RGBA. Free = translucent grey, bound = orange,
# bridge (couple holding two filaments) = red.
POINT_STATE_COLOR = {
    "free":   (0.55, 0.55, 0.60, 0.55),
    "bound":  (0.95, 0.55, 0.10, 1.0),
    "bridge": (0.90, 0.10, 0.10, 1.0),
    # an ATTACHED single exports two points: its anchor (e.g. on the membrane)
    # and the hand gripping the filament.  Different colours so the link stretch
    # between them is readable at a glance.
    "anchor": (0.20, 0.70, 0.95, 1.0),
    "hand":   (0.95, 0.35, 0.10, 1.0),
}

# Membrane appearance.  Semi-transparent by default: a solid shell hides every
# filament inside it, which is almost never what you want to look at.
MEMBRANE_COLOR = (0.60, 0.68, 0.80, 0.25)
# ===========================================================================


def _sorted_frames(directory, pattern):
    files = sorted(glob.glob(os.path.join(directory, pattern)))
    return files


def load_membrane_series(directory):
    """Return {label: [path per frame]} for every membrane sequence found.

    Two naming conventions coexist:
      membrane_%04d.obj      one membrane   (coupled_driver)
      membrane<K>_%04d.obj   K-th of N      (multi_membrane_driver)

    A single glob cannot cover both: "membrane_*.obj" does NOT match
    "membrane0_0000.obj", because the character after "membrane" is a digit, not
    the underscore. Multi-membrane runs were silently invisible in the viewer
    for exactly that reason -- the files were written, nothing was imported.
    """
    series = {}
    single = _sorted_frames(directory, "membrane_*.obj")
    if single:
        series["membrane"] = single
    for path in _sorted_frames(directory, "membrane[0-9]*_*.obj"):
        label = os.path.basename(path).rsplit("_", 1)[0]   # "membrane0"
        series.setdefault(label, []).append(path)
    for k in series:
        series[k].sort()
    if not series:
        print("[import_trajectory] no membrane OBJs in", directory)
    return series


def load_membrane_frames(directory):
    """Back-compat single-sequence accessor."""
    files = _sorted_frames(directory, "membrane_*.obj")
    if not files:
        print("[import_trajectory] no membrane OBJs in", directory)
    return files


def parse_fiber_file(path):
    """Parse one fibers_%04d.txt into a list of (radius, [(x,y,z), ...])."""
    fibers = []
    cur_pts = None
    cur_radius = FIBER_BEVEL_FALLBACK
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("fiber"):
                # flush previous
                if cur_pts:
                    fibers.append((cur_radius, cur_pts))
                # "fiber <id> <radius> <nPoints>"
                tok = line.split()
                cur_radius = float(tok[2]) if len(tok) > 2 else FIBER_BEVEL_FALLBACK
                cur_pts = []
            else:
                x, y, z = (float(v) for v in line.split()[:3])
                if cur_pts is not None:
                    cur_pts.append((x, y, z))
        if cur_pts:
            fibers.append((cur_radius, cur_pts))
    return fibers


def parse_point_file(path):
    """Parse one points_%04d.txt into a list of (class_name, state, [(x,y,z)...]).

    Block format written by coupled_driver:
        class <kind> <name> <state> <count>
        x y z
        ...
    """
    blocks = []
    cur = None
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith("class"):
                tok = line.split()
                # "class <kind> <name> <state> <count>"
                kind, name, state = tok[1], tok[2], tok[3]
                cur = ("%s %s %s" % (kind, name, state), name, state, [])
                blocks.append(cur)
            elif cur is not None:
                x, y, z = (float(v) for v in line.split()[:3])
                cur[3].append((x, y, z))
    return [(full, name, state, pts) for (full, name, state, pts) in blocks]


def _point_material(state):
    """One shared material per binding state, created on first use."""
    key = "cytosim_point_%s" % state
    mat = bpy.data.materials.get(key)
    if mat is None:
        mat = bpy.data.materials.new(key)
        mat.use_nodes = True
        rgba = POINT_STATE_COLOR.get(state, (0.7, 0.7, 0.7, 1.0))
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf is not None:
            bsdf.inputs["Base Color"].default_value = rgba
            if rgba[3] < 1.0:
                # older Blender spells this "Alpha"; guard in case it moves
                if "Alpha" in bsdf.inputs:
                    bsdf.inputs["Alpha"].default_value = rgba[3]
                mat.blend_method = 'BLEND'
        mat.diffuse_color = rgba          # viewport solid shading
    return mat


def _dot_mesh(radius, state):
    """Get/create the shared sphere MESH that point classes instance.

    Built with bmesh rather than bpy.ops.mesh.primitive_ico_sphere_add: the
    operator needs a scene context and leaves a stray object behind, which is
    awkward inside a frame-change handler.  The mesh datablock is cached and
    shared by every class with the same radius+state, so rebuilding a frame
    allocates objects but never geometry.
    """
    key = "cytosim_dot_%s_%.5g" % (state, radius)
    mesh = bpy.data.meshes.get(key)
    if mesh is not None:
        return mesh
    mesh = bpy.data.meshes.new(key)
    bm = bmesh.new()
    try:
        bmesh.ops.create_icosphere(bm, subdivisions=1, radius=radius)
    except TypeError:
        # Blender <3.0 spells the size argument "diameter"
        bmesh.ops.create_icosphere(bm, subdivisions=1, diameter=radius)
    bm.to_mesh(mesh)
    bm.free()
    mesh.materials.append(_point_material(state))
    return mesh


def build_point_clouds(blocks, coll, scale):
    """One vertex-instanced mesh per class.

    Instancing rather than one sphere per single: a run places hundreds of
    singles and couples, and hundreds of real objects per frame makes scrubbing
    the timeline unusable.  This way each class costs one object regardless of
    how many points it holds.
    """
    for full, name, state, pts in blocks:
        if not pts:
            continue
        radius = POINT_RADIUS_OVERRIDE.get(full,
                 POINT_RADIUS_OVERRIDE.get(name, POINT_RADIUS)) * scale
        mesh = bpy.data.meshes.new(full)
        mesh.from_pydata([(x * scale, y * scale, z * scale) for x, y, z in pts],
                         [], [])
        mesh.update()
        obj = bpy.data.objects.new(full, mesh)
        coll.objects.link(obj)
        # instance a sphere on every vertex.  The child must be its own object
        # (an object has exactly one parent), but the MESH behind it is shared.
        dot = bpy.data.objects.new(full + " [dot]", _dot_mesh(radius, state))
        coll.objects.link(dot)
        dot.parent = obj
        obj.instance_type = 'VERTS'


# distinct hues for membrane0, membrane1, ... so a multi-organelle scene is
# readable at a glance; index parsed from the label's trailing digits
_MEMBRANE_COLORS = [
    (0.45, 0.62, 0.90, 1.0),
    (0.90, 0.62, 0.40, 1.0),
    (0.55, 0.85, 0.55, 1.0),
    (0.85, 0.55, 0.80, 1.0),
]


def _membrane_material(label):
    key = "cytosim_membrane_%s" % label
    mat = bpy.data.materials.get(key)
    if mat is None:
        digits = "".join(c for c in label if c.isdigit())
        if digits:
            rgba = _MEMBRANE_COLORS[int(digits) % len(_MEMBRANE_COLORS)]
        else:
            rgba = MEMBRANE_COLOR          # the single-membrane case
        mat = bpy.data.materials.new(key)
        mat.use_nodes = True
        bsdf = mat.node_tree.nodes.get("Principled BSDF")
        if bsdf is not None:
            bsdf.inputs["Base Color"].default_value = rgba
            if rgba[3] < 1.0 and "Alpha" in bsdf.inputs:
                bsdf.inputs["Alpha"].default_value = rgba[3]
        if rgba[3] < 1.0:
            mat.blend_method = 'BLEND'
            # show the far side of the shell too, so the membrane reads as a
            # closed surface rather than a flat disc
            mat.use_backface_culling = False
        mat.diffuse_color = rgba
    return mat


def _clear_collection(coll):
    for obj in list(coll.objects):
        bpy.data.objects.remove(obj, do_unlink=True)


def build_fiber_curves(fibers, coll, scale):
    """Create one poly curve object per fiber in `fibers` into collection."""
    for i, (radius, pts) in enumerate(fibers):
        if len(pts) < 2:
            continue
        cu = bpy.data.curves.new("fiber_%d" % i, type='CURVE')
        cu.dimensions = '3D'
        spline = cu.splines.new('POLY')
        spline.points.add(len(pts) - 1)          # one point already exists
        for j, (x, y, z) in enumerate(pts):
            spline.points[j].co = (x * scale, y * scale, z * scale, 1.0)
        bevel = (radius * scale) if FIBER_BEVEL_FROM_RADIUS and radius > 0 \
            else FIBER_BEVEL_FALLBACK * scale
        cu.bevel_depth = bevel
        obj = bpy.data.objects.new("fiber_%d" % i, cu)
        coll.objects.link(obj)


def _get_or_make_collection(name):
    if name in bpy.data.collections:
        return bpy.data.collections[name]
    coll = bpy.data.collections.new(name)
    bpy.context.scene.collection.children.link(coll)
    return coll


def import_membrane_obj(path, scale):
    """Import one OBJ and return the created mesh object (or None).

    The axis arguments are REQUIRED, not cosmetic.  Blender's OBJ importer
    assumes the OBJ convention (Y up, -Z forward) and rotates the mesh into
    Blender's Z-up world.  Our fiber curves are built from raw simulation
    coordinates with no such conversion, so leaving the default on rotates the
    membrane 90 degrees about X relative to the filaments -- the protrusions
    then appear nowhere near the filament that caused them.

    Forcing forward=Y / up=Z makes the import an identity transform, so both
    objects share the one coordinate frame promised in COORDINATES.md.
    """
    before = set(bpy.data.objects)
    # Blender 4.x uses wm.obj_import; older uses import_scene.obj
    try:
        try:
            bpy.ops.wm.obj_import(filepath=path, forward_axis='Y', up_axis='Z')
        except TypeError:
            # older 4.x spelling of the same options
            bpy.ops.wm.obj_import(filepath=path, forward_axis='POSITIVE_Y',
                                  up_axis='POSITIVE_Z')
    except AttributeError:
        bpy.ops.import_scene.obj(filepath=path, axis_forward='Y', axis_up='Z')
    new = list(set(bpy.data.objects) - before)
    if not new:
        return None
    obj = new[0]
    obj.scale = (scale, scale, scale)
    return obj


# --- animation state, populated by setup() ---------------------------------
_STATE = {
    "membrane_series": {},
    "membrane_files": [],
    "fiber_files": [],
    "point_files": [],
    "membrane_coll": None,
    "fiber_coll": None,
    "point_coll": None,
}


def show_frame(frame_index):
    """Rebuild the scene for the given 0-based export frame index."""
    mem_series = _STATE["membrane_series"]
    fib_files = _STATE["fiber_files"]
    pt_files = _STATE["point_files"]
    mem_coll = _STATE["membrane_coll"]
    fib_coll = _STATE["fiber_coll"]
    pt_coll = _STATE["point_coll"]

    # membranes: clear, then import this frame of every sequence.  One
    # multi-membrane run has several, and they must all land in the same frame.
    if mem_coll is not None and mem_series:
        _clear_collection(mem_coll)
        for label in sorted(mem_series):
            files = mem_series[label]
            if not (0 <= frame_index < len(files)):
                continue
            obj = import_membrane_obj(files[frame_index], WORLD_SCALE)
            if obj is None:
                continue
            obj.name = label
            for c in obj.users_collection:
                c.objects.unlink(obj)
            mem_coll.objects.link(obj)
            # ALWAYS assign a material.  Previously this only happened when
            # there was more than one membrane, so a single-membrane run got an
            # untextured opaque mesh that simply hid the filaments inside it.
            obj.data.materials.append(_membrane_material(label))

    # fibers: clear then rebuild the curves for this frame
    if fib_coll is not None and 0 <= frame_index < len(fib_files):
        _clear_collection(fib_coll)
        fibers = parse_fiber_file(fib_files[frame_index])
        build_fiber_curves(fibers, fib_coll, WORLD_SCALE)

    # singles/couples: clear then rebuild the instanced point clouds
    if pt_coll is not None and 0 <= frame_index < len(pt_files):
        _clear_collection(pt_coll)
        blocks = parse_point_file(pt_files[frame_index])
        build_point_clouds(blocks, pt_coll, WORLD_SCALE)


def _frame_change_handler(scene):
    # map Blender timeline frame (1-based) to export index (0-based)
    idx = scene.frame_current - 1
    show_frame(idx)


def setup():
    _STATE["membrane_series"] = load_membrane_series(MEMBRANE_DIR)
    _STATE["membrane_files"] = next(
        iter(sorted(_STATE["membrane_series"].items())), ("", []))[1]
    _STATE["fiber_files"] = _sorted_frames(FIBER_DIR, "fibers_*.txt")
    _STATE["point_files"] = _sorted_frames(POINT_DIR, "points_*.txt")

    _STATE["membrane_coll"] = _get_or_make_collection("Membrane")
    _STATE["fiber_coll"] = _get_or_make_collection("Fibers")
    _STATE["point_coll"] = _get_or_make_collection("Points")

    nframes = max(len(_STATE["membrane_files"]), len(_STATE["fiber_files"]),
                  len(_STATE["point_files"]))
    if nframes == 0:
        print("[import_trajectory] nothing to import -- check the paths in CONFIG")
        return

    scene = bpy.context.scene
    scene.frame_start = 1
    scene.frame_end = nframes

    # register the handler (remove any stale copy first)
    handlers = bpy.app.handlers.frame_change_post
    handlers[:] = [h for h in handlers if h.__name__ != "_frame_change_handler"]
    handlers.append(_frame_change_handler)

    scene.frame_set(1)
    print("[import_trajectory] ready: %d membrane, %d fiber, %d point frame(s)"
          % (len(_STATE["membrane_files"]), len(_STATE["fiber_files"]),
             len(_STATE["point_files"])))
    _check_registration()


def _check_registration():
    """Warn if the imported membrane does not match the file's own coordinates.

    Catches the classic OBJ axis-conversion trap: if Blender rotated the mesh on
    import, the membrane and the fiber curves end up in different frames and the
    protrusions appear detached from the filaments that made them.
    """
    files = _STATE["membrane_files"]
    coll = _STATE["membrane_coll"]
    if not files or coll is None or not coll.objects:
        return

    # first vertex as written in the OBJ
    want = None
    with open(files[0]) as fh:
        for line in fh:
            if line.startswith("v "):
                want = tuple(float(t) for t in line.split()[1:4])
                break
    if want is None:
        return

    obj = coll.objects[0]
    got = tuple(obj.matrix_world @ obj.data.vertices[0].co)
    err = max(abs(a * WORLD_SCALE - b) for a, b in zip(want, got))
    if err > 1e-4 * max(1.0, max(abs(v) for v in want)):
        print("[import_trajectory] WARNING: membrane was transformed on import")
        print("    file says   %s" % (want,))
        print("    Blender has %s" % (got,))
        print("    -> the membrane and fibers are in different frames;")
        print("       check the forward_axis/up_axis arguments in "
              "import_membrane_obj()")
    else:
        print("[import_trajectory] registration OK "
              "(membrane matches file coordinates)")


if __name__ == "__main__":
    setup()
