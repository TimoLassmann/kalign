const std = @import("std");

const kalignPackageVersion = "3.5.1";

const targets: []const std.Target.Query = &.{
    .{ .cpu_arch = .aarch64, .os_tag = .macos },
    .{ .cpu_arch = .aarch64, .os_tag = .linux },
    .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .gnu },
    .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .musl },
};

const cflags = [_][]const u8{
    "-DKALIGN_PACKAGE_VERSION=\"" ++ kalignPackageVersion ++ "\"",
    "-DKALIGN_PACKAGE_NAME=\"kalign\"",
    "-DKALIGN_ALN_SERIAL_THRESHOLD=250",
    "-DKALIGN_KMEANS_UPGMA_THRESHOLD=50",
};

pub fn build(b: *std.Build) !void {
    const optimize = b.standardOptimizeOption(.{});

    for (targets) |t| {
        // --- Static library module ---
        const lib_mod = b.createModule(.{
            .target = b.resolveTargetQuery(t),
            .optimize = optimize,
            .link_libc = true,
        });
        lib_mod.addIncludePath(b.path("lib/src"));
        lib_mod.addIncludePath(b.path("lib/include"));
        lib_mod.addCSourceFiles(.{ .files = &kalign_lib_sources, .flags = &cflags });

        const lib = b.addLibrary(.{
            .name = "tldevel",
            .linkage = .static,
            .root_module = lib_mod,
        });
        b.installArtifact(lib);

        // --- Executable module ---
        const bin_mod = b.createModule(.{
            .target = b.resolveTargetQuery(t),
            .optimize = optimize,
            .link_libc = true,
        });
        bin_mod.addIncludePath(b.path("lib/src"));
        bin_mod.addIncludePath(b.path("lib/include"));
        bin_mod.addCSourceFiles(.{ .files = &kalign_sources, .flags = &cflags });
        bin_mod.linkLibrary(lib);

        const kalign_bin = b.addExecutable(.{
            .name = "kalign",
            .root_module = bin_mod,
        });
        b.installArtifact(kalign_bin);

        // Install into a per-target subdirectory (e.g. zig-out/x86_64-linux-gnu/kalign)
        const target_output = b.addInstallArtifact(kalign_bin, .{
            .dest_dir = .{
                .override = .{
                    .custom = try t.zigTriple(b.allocator),
                },
            },
        });

        b.getInstallStep().dependOn(&target_output.step);
    }
}

const kalign_lib_sources = [_][]const u8{
    "lib/src/test.c",
    "lib/src/tldevel.c",
    "lib/src/tlmisc.c",
    "lib/src/tlrng.c",
    "lib/src/esl_stopwatch.c",
    "lib/src/msa_alloc.c",
    "lib/src/msa_op.c",
    "lib/src/msa_io.c",
    "lib/src/msa_misc.c",
    "lib/src/msa_check.c",
    "lib/src/msa_cmp.c",
    "lib/src/msa_sort.c",
    "lib/src/alphabet.c",
    "lib/src/task.c",
    "lib/src/bisectingKmeans.c",
    "lib/src/sequence_distance.c",
    "lib/src/bpm.c",
    "lib/src/euclidean_dist.c",
    "lib/src/pick_anchor.c",
    "lib/src/aln_wrap.c",
    "lib/src/aln_apair_dist.c",
    "lib/src/aln_add.c",
    "lib/src/aln_param.c",
    "lib/src/aln_run.c",
    "lib/src/aln_mem.c",
    "lib/src/aln_setup.c",
    "lib/src/aln_controller.c",
    "lib/src/aln_seqseq.c",
    "lib/src/aln_seqprofile.c",
    "lib/src/aln_profileprofile.c",
    "lib/src/aln_refine.c",
    "lib/src/sp_score.c",
    "lib/src/weave_alignment.c",
    "lib/src/poar.c",
    "lib/src/consensus_msa.c",
    "lib/src/anchor_consistency.c",
    "lib/src/msa_consistency.c",
    "lib/src/ensemble.c",
};

const kalign_sources = [_][]const u8{
    "src/run_kalign.c",
    "src/parameters.c",
};
