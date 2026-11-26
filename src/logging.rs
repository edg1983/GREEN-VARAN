use anyhow::{Context, Result};
use chrono::Utc;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::time::Instant;

pub struct DualLogger {
    file: Option<Mutex<BufWriter<File>>>,
}

impl DualLogger {
    pub fn console_only() -> Self {
        Self { file: None }
    }

    pub fn new(log_opt: Option<&str>, out_path: &str) -> Result<Self> {
        let log_path: PathBuf = match log_opt {
            Some(p) if !p.is_empty() => PathBuf::from(p),
            _ => {
                let ts = Utc::now().format("%Y%m%d_%H%M%S");
                let parent = Path::new(out_path)
                    .parent()
                    .filter(|p| !p.as_os_str().is_empty())
                    .unwrap_or_else(|| Path::new("."));
                parent.join(format!("greenvaran_{ts}.log"))
            }
        };
        let file = File::create(&log_path)
            .with_context(|| format!("Failed to create log file: {}", log_path.display()))?;
        Ok(Self {
            file: Some(Mutex::new(BufWriter::new(file))),
        })
    }

    fn format_line(level: &str, msg: &str) -> String {
        let now = Utc::now().format("%Y-%m-%dT%H:%M:%S");
        format!("[{now}] - {level}: {msg}\n")
    }

    fn write_file(&self, line: &str) {
        if let Some(ref f) = self.file {
            if let Ok(mut guard) = f.lock() {
                let _ = guard.write_all(line.as_bytes());
                let _ = guard.flush();
            }
        }
    }

    pub fn debug(&self, msg: &str) {
        let line = Self::format_line("DEBUG", msg);
        self.write_file(&line);
    }

    pub fn info(&self, msg: &str) {
        let line = Self::format_line("INFO", msg);
        print!("{line}");
        self.write_file(&line);
    }

    pub fn warn(&self, msg: &str) {
        let line = Self::format_line("WARN", msg);
        print!("{line}");
        self.write_file(&line);
    }

    pub fn error(&self, msg: &str) {
        let line = Self::format_line("ERROR", msg);
        eprint!("{line}");
        self.write_file(&line);
    }

    pub fn fatal(&self, msg: &str) {
        let line = Self::format_line("FATAL", msg);
        eprint!("{line}");
        self.write_file(&line);
    }
}

pub fn format_elapsed(t0: &mut Instant) -> String {
    let elapsed = t0.elapsed().as_secs_f64();
    let mins = (elapsed / 60.0).floor();
    let secs = elapsed - mins * 60.0;
    *t0 = Instant::now();
    format!("{mins:.1} min {secs:.2} sec")
}
