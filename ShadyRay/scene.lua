--Sample scene description file for the ShadyRay raytracer.

--Custom functions exported by ShadyRay: render(render_ops, materials) (don't look for this in the Lua libraries!

--You can use all Lua functionality. IMPORTANT: the Lua script's last action should be to call render().
--Configuration is done by deepcopying Lua tables containing information in a specific structure.
--For convenience several default tables have been provided. deepcopy these and change what you need.

--Use this function to copy tables. IMPORTANT: do not simply assign tables! Lua tables are based on references.
function deepcopy(orig)
    local orig_type = type(orig)
    local copy

    if orig_type == 'table' then
        copy = {}

        for orig_key, orig_value in next, orig, nil do
            copy[deepcopy(orig_key)] = deepcopy(orig_value)
        end

        setmetatable(copy, deepcopy(getmetatable(orig)))

    else
        copy = orig

    end

    return copy
end

--This is a standard material. IMPORTANT: The mesh name must match that exported by the 3D application.
DEFAULT_MATERIAL = {
    mesh_name = "Mesh",
    diffuse_color = {
        r = 1.0,
        g = 1.0,
        b = 1.0,
    },

    specular_color = {
        r = 0.0,
        g = 0.0,
        b = 0.0,
    },

    transmission_color = {
        r = 0.0,
        g = 0.0,
        b = 0.0,

        ior = 1.1
    },

    emission_color = {
        r = 0.0,
        g = 0.0,
        b = 0.0
    },

    alpha = 0.01,

    color_map = {
        type = "Constant",
        path = "",
    },

    normal_map = {
        type = "Constant",
        path = "",
    },

    mapping = "PlaneY",

    is_lamp = 0
}

--This is a standard material table.
DEFAULT_MATERIALS = {
    num_materials = 0,
    material_table = {}
}


--This is a standard scene description. Change the path to your asset file.
DEFAULT_OPTIONS = {
    scene = {
        path = "/home/walksbynight/Assets/transparent.dae",
        width = 1024,
        height = 768,
        scale = 1.0,
    },

    camera = {
        focal_length = 4.0,

        scale = {
            x = 0.003,
            y = 0.003
        },

        origin = {
            x = 2.0,
            y = 2.0,
            z = 5
        },

        rotation = {
            x = 0,
            y = 0.5,
            z = 0
        }
    },

    sampler = {
        multiplier = 2.0,
        num_samples = 8
    },

    filter = {
        b = 0.3,
        c = 0.3,
        radius = 2.0
    },

    integrator = {
        x_blocks = 16,
        y_blocks = 16,
        num_bounces = 8
    },

    window_name = "Lua Test"
}


--Sample usage

--First copy the material table and set the number of materials.
materials = deepcopy(DEFAULT_MATERIALS)
materials.num_materials = 3

--Then insert materials in the material table. IMPORTANT: must use string keys, not numeral. Starts at "1".
materials.material_table["1"] = deepcopy(DEFAULT_MATERIAL)
materials.material_table["1"].mesh_name = "Icosphere.001"
materials.material_table["1"].emission_color.r = 1.0
materials.material_table["1"].emission_color.g = 1.0
materials.material_table["1"].emission_color.b = 1.0


materials.material_table["2"] = deepcopy(DEFAULT_MATERIAL)
materials.material_table["2"].mesh_name = "Untitled.017"
materials.material_table["2"].diffuse_color = {r = 0, g = 0, b = 0 }
materials.material_table["2"].alpha = 0.1
materials.material_table["2"].transmission_color = {r = 1, g = 1, b = 1, ior = 1.7}
materials.material_table["2"].is_lamp = 1


materials.material_table["3"] = deepcopy(DEFAULT_MATERIAL)
materials.material_table["3"].mesh_name = "Plane.001"

--Then copy the scene description
render_options = deepcopy(DEFAULT_OPTIONS)

--Make changes to scene description
render_options.sampler.num_samples = 1024
render_options.integrator.num_bounces = 8


--Call the render function
render(render_options, materials)