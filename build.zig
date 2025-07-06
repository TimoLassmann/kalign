const std = @import("std");
const builtin = @import("builtin");
const ArrayList = std.ArrayList;

const kalignPackageVersion = "3.4.1";

const targets: []const std.Target.Query = &.{
    .{ .cpu_arch = .aarch64, .os_tag = .macos },
    .{ .cpu_arch = .aarch64, .os_tag = .linux },
    .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .gnu },
    .{ .cpu_arch = .x86_64, .os_tag = .linux, .abi = .musl },
    // .{ .cpu_arch = .x86_64, .os_tag = .windows },
};

const cflags = [_][]const u8{
    "-DKALIGN_PACKAGE_VERSION=\"3.4.1\"",
    "-DKALIGN_PACKAGE_NAME=\"kalign\"",
    "-DKALIGN_ALN_SERIAL_THRESHOLD=250",
    "-DKALIGN_KMEANS_UPGMA_THRESHOLD=50",
};

pub fn build(b: *std.Build) !void {
    // const target = b.standardTargetOptions(.{});
    const optimize = b.standardOptimizeOption(.{});

    for (targets) |t| {
        const lib = b.addStaticLibrary(.{
            .name = "tldevel",
            .target = b.resolveTargetQuery(t),
            .optimize = optimize,
        });
        lib.linkLibC();
        lib.addIncludePath(.{ .path = "./lib/src" });
        lib.addCSourceFiles(.{ .files = &kalign_lib_sources, .flags = &cflags });
        lib.addIncludePath(.{ .path = "./lib/include" });
        b.installArtifact(lib);

        const kalign_bin = b.addExecutable(.{
            .name = "kalign",
            .target = b.resolveTargetQuery(t),
            .optimize = optimize,
        });

        lib.addCSourceFiles(.{ .files = &kalign_sources, .flags = &cflags });
        kalign_bin.addIncludePath(.{ .path = "./lib/src" });
        kalign_bin.addIncludePath(.{ .path = "./lib/include" });
        kalign_bin.linkLibrary(lib);
        b.installArtifact(kalign_bin);

        const target_output = b.addInstallArtifact(kalign_bin, .{
            .dest_dir = .{
                .override = .{
                    .custom = try t.zigTriple(b.allocator),
                },
            },
        });

        b.getInstallStep().dependOn(&target_output.step);
    }

    // const exe = b.addExecutable(.{
    //     .name = "zig_test",
    //     .target = target,
    //     .optimize = optimize,
    // });
    // exe.addCSourceFile(.{ .file = .{ .path = "./tests/zig_test.c" }, .flags = &[_][]const u8{"-std=c99"} });

    // // exe.addCSourceFiles(.{ .files = &kalign_lib_sources, .flags = &[_][]const u8{"-std=c99"} });
    // // exe.addCSourceFile(&.{"./tests/zig_test.c"}, cflags.items);
    // // exe.addCSourceFile("./lib/src/strnlen_compat.c", cflags.items);
    // exe.addIncludePath(.{ .path = "./lib/src" });
    // // exe.linkLibrary(lib);
    // exe.linkLibC();
    // b.installArtifact(exe);
}

const kalign_lib_sources = [_][]const u8{
    // "lib/src/strnlen_compat.c",
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
    "lib/src/aln_param.c",
    "lib/src/aln_run.c",
    "lib/src/aln_mem.c",
    "lib/src/aln_setup.c",
    "lib/src/aln_controller.c",
    "lib/src/aln_seqseq.c",
    "lib/src/aln_seqprofile.c",
    "lib/src/aln_profileprofile.c",
    "lib/src/weave_alignment.c",
};

const kalign_sources = [_][]const u8{
    "./src/run_kalign.c",
    "./src/parameters.c",
};
