using Pkg

println("[NetInfer] Julia version info before install:")
try
	import InteractiveUtils

	# Hardened install with retries and fallbacks for environments with curl issues
	try
		# Prefer using git instead of libcurl-based downloader when environments are problematic
		ENV["JULIA_PKG_USE_CLI_GIT"] = get(ENV, "JULIA_PKG_USE_CLI_GIT", "true")
		# Disable package servers if they cause issues; fallback to direct Git
		ENV["JULIA_PKG_SERVER"] = get(ENV, "JULIA_PKG_SERVER", "")

		println("[NetInfer] Julia version info:")
		InteractiveUtils.versioninfo()
		println("[NetInfer] Installing FlashWeave.jl …")

		success = false
		max_tries = 3
		for attempt in 1:max_tries
			try
				Pkg.add("FlashWeave")
				Pkg.precompile()
				success = true
				break
			catch e
				@warn "FlashWeave install attempt $attempt failed: $e"
				sleep(5)
			end
		end

		if !success
			error("Failed to install FlashWeave after $(max_tries) attempts.")
		end

		println("[NetInfer] FlashWeave installation completed.")
	catch e
		@error "FlashWeave installation error: $e"
		rethrow()
	end
catch e
	@warn "versioninfo() failed" exception=(e, catch_backtrace())
end

# Harden downloads in restricted environments:
# - Prefer system curl if bundled one misbehaves
# - Optionally disable Pkg servers to fetch directly (can help behind proxies/mitm)
ENV["DOWNLOADS_USESYSTEMCURL"] = get(ENV, "DOWNLOADS_USESYSTEMCURL", "1")
ENV["JULIA_PKG_SERVER"] = get(ENV, "JULIA_PKG_SERVER", "")  # empty disables pkg server

try
	println("[NetInfer] Adding FlashWeave (first attempt)...")
	Pkg.add("FlashWeave")
catch e
	@warn "FlashWeave add failed; retrying with registry update and no Pkg server" exception=(e, catch_backtrace())
	try
		Pkg.Registry.update()
	catch e2
		@warn "Registry update failed" exception=(e2, catch_backtrace())
	end
	# Force disable Pkg server on retry
	ENV["JULIA_PKG_SERVER"] = ""
	Pkg.add("FlashWeave")
end

println("[NetInfer] FlashWeave installation complete.")