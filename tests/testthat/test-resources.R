# THREADS=0 resolves from the cgroup limits (fake cgroup trees), not the host's CPU count.

fake_cgroup <- function(files) {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  for (n in names(files)) {
    f <- file.path(root, n)
    dir.create(dirname(f), recursive = TRUE, showWarnings = FALSE)
    writeLines(files[[n]], f)
  }
  root
}

test_that("cgroup v2: the CPU quota, capped by the memory limit", {
  root <- fake_cgroup(list(cpu.max = "400000 100000", memory.max = as.character(3 * 2^30)))
  expect_equal(cgroup_cpus(root), 4)
  expect_equal(cgroup_memory_bytes(root), 3 * 2^30)
  w <- resolve_workers(0, 550, 500, root, affinity = 15L, machine = 64 * 2^30)
  expect_equal(w$n, 4L)
  expect_match(w$why, "cgroup CPU quota 4")
  w <- resolve_workers(0, 1000, 500, root, affinity = 15L, machine = 64 * 2^30)
  expect_equal(w$n, 2L)                    # (3072 - 500) %/% 1000
  expect_match(w$why, "capped by memory")
})

test_that("no limits: CPU affinity and the machine's memory", {
  root <- fake_cgroup(list(cpu.max = "max 100000", memory.max = "max"))
  expect_null(cgroup_cpus(root))
  expect_null(cgroup_memory_bytes(root))
  w <- resolve_workers(0, 550, 500, root, affinity = 15L, machine = 64 * 2^30)
  expect_equal(w$n, 15L)
  expect_match(w$why, "CPU affinity")
})

test_that("cgroup v1: fractional quota, unlimited memory, -1 quota", {
  root <- fake_cgroup(list("cpu/cpu.cfs_quota_us" = "150000", "cpu/cpu.cfs_period_us" = "100000",
                           "memory/memory.limit_in_bytes" = "9223372036854771712"))
  expect_equal(cgroup_cpus(root), 1.5)
  expect_null(cgroup_memory_bytes(root))
  expect_equal(resolve_workers(0, 550, 500, root, affinity = 15L, machine = 64 * 2^30)$n, 1L)
  root2 <- fake_cgroup(list("cpu/cpu.cfs_quota_us" = "-1", "cpu/cpu.cfs_period_us" = "100000"))
  expect_null(cgroup_cpus(root2))
})

test_that("an explicit THREADS is kept, with a note when memory looks short", {
  root <- fake_cgroup(list(cpu.max = "200000 100000", memory.max = as.character(2 * 2^30)))
  w <- resolve_workers(8, 550, 500, root, affinity = 15L, machine = 64 * 2^30)
  expect_equal(w$n, 8L)
  expect_match(w$why, "more than the 2")
  expect_equal(resolve_workers(2, 550, 500, root, affinity = 15L)$why, "THREADS=2")
})

test_that("the affinity mask is read from /proc/self/status", {
  f <- withr::local_tempfile()
  writeLines(c("Name:\tR", "Cpus_allowed:\tff0f", "Cpus_allowed_list:\t0-3,8-11,14"), f)
  expect_equal(affinity_cpus(f), 9L)
  expect_null(affinity_cpus(file.path(tempdir(), "missing")))
})
