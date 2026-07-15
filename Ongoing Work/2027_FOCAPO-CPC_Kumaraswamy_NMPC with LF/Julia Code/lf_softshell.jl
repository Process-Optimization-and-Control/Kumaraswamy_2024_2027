function resolve_softsh_path()
    if haskey(ENV, "LEDAFLOW_SOFTSH")
        softsh_env = ENV["LEDAFLOW_SOFTSH"]
        isfile(softsh_env) && return softsh_env
        error("LEDAFLOW_SOFTSH is set but file was not found: $(softsh_env)")
    end

    candidates = Sys.iswindows() ?
        ["C:/Program Files/Kongsberg/LedaFlow Engineering v2.11.271.018/softsh.exe"] :
        ["/mnt/c/Program Files/Kongsberg/LedaFlow Engineering v2.11.271.018/softsh.exe"]

    for p in candidates
        isfile(p) && return p
    end

    error("Could not find LedaFlow softsh.exe. Set LEDAFLOW_SOFTSH to the full executable path.")
end