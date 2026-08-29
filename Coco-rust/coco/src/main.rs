//! CoCo command dispatch, ported from `main.cpp`.

use coco::commands;
use coco::info::Info;
use coco::options::{parse_options, CommandId, ParseOutcome};

struct CommandSpec {
    cmd: &'static str,
    id: CommandId,
    descript_short: &'static str,
    author: &'static str,
    usage: &'static str,
}

const COMMANDS: &[CommandSpec] = &[
    CommandSpec {
        cmd: "correction",
        id: CommandId::Correction,
        descript_short: "correct sequencing errors",
        author: "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        usage: " -1 <fasta|q> -2 <fasta|q> | --reads <fasta|q> [--counts <count.h5>] [options]",
    },
    CommandSpec {
        cmd: "filter",
        id: CommandId::Filter,
        descript_short: "filter chimeric reads",
        author: "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        usage: " -1 <fasta|q> -2 <fasta|q> | --reads <fasta|q> [--counts <count.h5>] [options]",
    },
    CommandSpec {
        cmd: "abundance",
        id: CommandId::Abundance,
        descript_short: "estimate abundance values",
        author: "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        usage: "  -1 <fasta|q> -2 <fasta|q> | --reads <fasta|q> [--counts <count.h5>] [options]",
    },
    CommandSpec {
        cmd: "profile",
        id: CommandId::Profile,
        descript_short: "print spaced k-mer count profiles (devtool)",
        author: "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        usage: "  -1 <fasta|q> -2 <fasta|q> | --reads <fasta|q> [--counts <count.h5>] [options]",
    },
    CommandSpec {
        cmd: "counts2flat",
        id: CommandId::Counts2Flat,
        descript_short: "print spaced k-mer lookup table (devtool)",
        author: "Annika Jochheim <annika.jochheim@mpinat.mpg.de>",
        usage: " --counts <count.h5>] [--outprefix <string>] [options]",
    },
];

fn print_usage() {
    let mut usage = String::new();
    usage.push_str(coco::TOOL_INTRODUCTION);
    usage.push_str("\n\n");
    usage.push_str(&format!("{} Version: {}\n", coco::TOOL_NAME, coco::VERSION));
    usage.push_str(&format!("© {}\n\n", coco::MAIN_AUTHOR));
    usage.push_str("available commands:\n");
    for c in COMMANDS {
        usage.push_str(&format!("  {}\t{}\n", c.cmd, c.descript_short));
    }
    eprintln!("{usage}");
}

fn main() {
    let argv: Vec<String> = std::env::args().collect();
    if argv.len() < 2 {
        print_usage();
        std::process::exit(1);
    }
    if argv[1] == "--help" || argv[1].starts_with("-h") {
        print_usage();
        std::process::exit(0);
    }

    let Some(spec) = COMMANDS.iter().find(|c| c.cmd == argv[1]) else {
        eprintln!("Invalid Command: {}", argv[1]);
        print_usage();
        std::process::exit(1);
    };

    // The C++ prints this from main(), before the command parses its options --
    // and before Info::setVerboseLevel runs, so it appears even at --verbose 0.
    let startup = Info::new(Info::INFO);
    startup.info(&format!("\n{} Version: {}\n", coco::TOOL_NAME, coco::VERSION));
    startup.info(&format!("Execute {} command: {}\n\n", coco::TOOL_NAME, spec.cmd));

    let rest: Vec<String> = argv[2..].to_vec();
    let opt = match parse_options(&rest, spec.id, spec.author, spec.usage, spec.cmd) {
        ParseOutcome::Ok(o) => o,
        ParseOutcome::HelpRequested(text) => {
            eprint!("{text}");
            std::process::exit(0);
        }
        ParseOutcome::Error(e) => {
            eprint!("{e}");
            std::process::exit(1);
        }
    };

    let info = Info::new(opt.verbose);
    info.info(&opt.parameter_settings(spec.id));
    info.info("\n");

    if opt.threads > 1 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(opt.threads as usize)
            .build_global()
            .ok();
    }

    let res = match spec.id {
        CommandId::Correction => commands::correction::run(&opt, &info),
        CommandId::Filter => commands::filter::run(&opt, &info),
        CommandId::Abundance => commands::abundance::run(&opt, &info),
        CommandId::Profile => commands::profile::run(&opt, &info),
        CommandId::Counts2Flat => commands::counts2flat::run(&opt, &info),
    };

    match res {
        Ok(code) => std::process::exit(code),
        Err(e) => {
            info.error(&e);
            std::process::exit(1);
        }
    }
}
