mod calc;
mod cli;
mod error;

fn main() -> std::result::Result<(), ()> {
    match main_try() {
        Ok(()) => Ok(()),
        Err(e) => {
            eprintln!("{}", e);
            Err(())
        }
    }
}

fn main_try() -> Result<(), Box<dyn std::error::Error>> {
    Ok(())
}
