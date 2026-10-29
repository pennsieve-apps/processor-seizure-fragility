# How many CPUs and how much memory this container may use: its cgroup limits, not the host's.
#
# parallel::detectCores() reports the host (or VM) CPUs, which inside a CPU-limited container
# (Docker --cpus, an ECS/Fargate task) can be several times the allocation; forking that many
# workers under a memory limit gets the run OOM-killed without a word. THREADS=0 is therefore
# resolved from the cgroup CPU quota (v2 cpu.max, v1 cpu.cfs_quota_us / cpu.cfs_period_us), else
# the CPU affinity (nproc), else detectCores(), and capped so that the estimated per-worker
# memory fits the cgroup memory limit (v2 memory.max, v1 memory.limit_in_bytes), else the
# machine's memory. Same rules as processor/resources.py of the Python epilepsy iEEG apps.

CGROUP <- "/sys/fs/cgroup"

read_first <- function(path) {
  if (!file.exists(path)) return(NULL)
  v <- tryCatch(trimws(readLines(path, n = 1, warn = FALSE)), error = function(e) NULL)
  if (!length(v)) NULL else v
}

# CPUs allowed by the cgroup CPU quota, or NULL when there is no quota.
cgroup_cpus <- function(root = CGROUP) {
  v2 <- read_first(file.path(root, "cpu.max"))
  if (!is.null(v2)) {
    parts <- strsplit(v2, " ")[[1]]
    if (parts[1] == "max") return(NULL)
    q <- suppressWarnings(as.numeric(parts[1]))
    p <- suppressWarnings(as.numeric(if (length(parts) > 1) parts[2] else 100000))
    return(if (is.na(q) || is.na(p) || p <= 0) NULL else q / p)
  }
  for (d in file.path(root, c("cpu", "cpu,cpuacct", "."))) {
    q <- read_first(file.path(d, "cpu.cfs_quota_us"))
    p <- read_first(file.path(d, "cpu.cfs_period_us"))
    if (!is.null(q) && !is.null(p)) {
      q <- suppressWarnings(as.numeric(q)); p <- suppressWarnings(as.numeric(p))
      return(if (is.na(q) || is.na(p) || q <= 0 || p <= 0) NULL else q / p)
    }
  }
  NULL
}

# The cgroup memory limit in bytes, or NULL when there is none.
cgroup_memory_bytes <- function(root = CGROUP) {
  for (f in file.path(root, c("memory.max", "memory/memory.limit_in_bytes", "memory.limit_in_bytes"))) {
    v <- read_first(f)
    if (is.null(v)) next
    if (v == "max") return(NULL)
    n <- suppressWarnings(as.numeric(v))
    return(if (is.na(n) || n <= 0 || n >= 2^60) NULL else n)
  }
  NULL
}

machine_memory_bytes <- function() {
  m <- read_first("/proc/meminfo")   # "MemTotal:  16314896 kB"
  if (is.null(m)) return(NULL)
  n <- suppressWarnings(as.numeric(gsub("[^0-9]", "", m)))
  if (is.na(n)) NULL else n * 1024
}

# CPUs in this process's affinity mask, from /proc/self/status Cpus_allowed_list ("0-3,8").
# (Not nproc: it honours OMP_NUM_THREADS, which images set for BLAS.)
affinity_cpus <- function(status = "/proc/self/status") {
  if (!file.exists(status)) return(NULL)
  l <- grep("^Cpus_allowed_list:", readLines(status, warn = FALSE), value = TRUE)
  if (!length(l)) return(NULL)
  n <- 0L
  for (r in strsplit(trimws(sub("^[^:]*:", "", l[1])), ",")[[1]]) {
    b <- suppressWarnings(as.integer(strsplit(r, "-")[[1]]))
    if (anyNA(b) || !length(b)) return(NULL)
    n <- n + if (length(b) == 2) b[2] - b[1] + 1L else 1L
  }
  if (n < 1) NULL else n
}

# Whole CPUs this process may use (at least 1) and where the number comes from.
available_cpus <- function(root = CGROUP, affinity = affinity_cpus()) {
  quota <- cgroup_cpus(root)
  if (!is.null(quota)) {
    n <- max(1L, as.integer(floor(quota + 1e-9)))
    if (!is.null(affinity)) n <- min(n, affinity)
    return(list(n = n, why = sprintf("cgroup CPU quota %g", quota)))
  }
  if (!is.null(affinity)) return(list(n = affinity, why = "CPU affinity"))
  list(n = max(1L, parallel::detectCores(), na.rm = TRUE), why = "detectCores()")
}

# Number of parallel workers and why. requested > 0 is used as given (with a note when the
# memory estimate says it will not fit); 0 = one per available CPU, capped so that reserve_mb
# (the parent) plus per_worker_mb per worker fits the memory limit.
resolve_workers <- function(requested, per_worker_mb, reserve_mb, root = CGROUP,
                            affinity = affinity_cpus(), machine = machine_memory_bytes()) {
  limit <- cgroup_memory_bytes(root)
  src <- "cgroup memory limit"
  if (is.null(limit)) { limit <- machine; src <- "machine memory" }
  fit <- if (is.null(limit)) NULL else max(1L, as.integer((limit / 2^20 - reserve_mb) %/% per_worker_mb))
  note <- if (is.null(fit)) "memory limit unknown" else
    sprintf("%s %.0f MB fits %d at ~%g MB each + %g MB", src, limit / 2^20, fit, per_worker_mb, reserve_mb)
  if (requested > 0) {
    why <- sprintf("THREADS=%d", requested)
    if (!is.null(fit) && requested > fit)
      why <- sprintf("%s (more than the %d the memory estimate allows: %s)", why, fit, note)
    return(list(n = as.integer(requested), why = why))
  }
  cpu <- available_cpus(root, affinity)
  if (!is.null(fit) && fit < cpu$n)
    return(list(n = fit, why = sprintf("THREADS=0: %d CPU(s) (%s), capped by memory: %s", cpu$n, cpu$why, note)))
  list(n = cpu$n, why = sprintf("THREADS=0: %d CPU(s) (%s); %s", cpu$n, cpu$why, note))
}
