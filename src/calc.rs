//// #### A little calculator parser and evaluator for any Field
////
//// Note: This source file will use some of the more advanced Rust features
//// to create a generic calculator that we will reuse for each new field we construct.
////
//// Understanding this code is not required nor expected to follow the series. This builds
//// a tokenizer and recursive-descent parser and is very much out of scope for the series.
//// If you'd like to read it, you're more than welcome. But don't spend too much time if it's
//// confusing. It's not written for learning as much as the Field implementations.
////
//// But, as a tool it's invaluable to use for studying those Fields!
use crate::field::Field;
use std::io::{self, BufRead, Write};

#[derive(Debug)]
pub enum Error {
  ParseFailed(String),
  InvalidExpression,
  DivisionByZero,
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Operator {
  Add,
  Sub,
  Mul,
  Div,
}

const MAX_PREC: usize = 2;

impl Operator {
  fn prec(&self) -> usize {
    match self {
      Operator::Add => 2,
      Operator::Sub => 2,
      Operator::Mul => 1,
      Operator::Div => 1,
    }
  }
}

#[derive(Debug, Clone, Copy, PartialEq)]
enum Token<F: Field> {
  Number(F),
  Oper(Operator),
  LParen,
  RParen,
}

struct Lexer<'a, F: Field> {
  buf: &'a [u8],
  idx: usize,
  tok: Option<Token<F>>,
}

impl<'a, F: Field> Lexer<'a, F> {
  fn new(input: &'a str) -> Result<Self, Error> {
    let mut this = Self {
      buf: input.as_bytes(),
      idx: 0,
      tok: None,
    };

    this.token_next()?;
    Ok(this)
  }

  fn peek_char(&self) -> Option<char> {
    if self.idx == self.buf.len() {
      return None;
    }
    Some(char::from_u32(self.buf[self.idx] as u32).unwrap())
  }

  fn advance_char(&mut self) {
    self.idx += 1;
  }

  fn token_next(&mut self) -> Result<(), Error> {
    loop {
      let mut ch = match self.peek_char() {
        Some(ch) => ch,
        None => {
          // End of text
          self.tok = None;
          return Ok(());
        }
      };

      // skip white space
      if ch.is_whitespace() {
        self.advance_char();
        continue;
      }

      // parse num?
      if ch.is_ascii_digit() {
        let start = self.idx;
        while ch.is_ascii_digit() {
          self.advance_char();
          ch = match self.peek_char() {
            Some(ch) => ch,
            None => break,
          };
        }
        let end = self.idx;

        let num_str = std::str::from_utf8(&self.buf[start..end]).unwrap();
        let num: F = num_str.parse().map_err(|s| Error::ParseFailed(s))?;

        self.tok = Some(Token::Number(num));
        return Ok(());
      }

      // parse various punctuation
      let tok = match ch {
        '+' => Token::Oper(Operator::Add),
        '-' => Token::Oper(Operator::Sub),
        '*' => Token::Oper(Operator::Mul),
        '/' => Token::Oper(Operator::Div),
        '(' => Token::LParen,
        ')' => Token::RParen,
        _ => return Err(Error::ParseFailed(format!("Invalid punctuation: '{}'", ch))),
      };

      self.advance_char();

      self.tok = Some(tok);
      return Ok(());
    }
  }

  fn token_peek(&self) -> Option<Token<F>> {
    self.tok
  }

  fn expect(&mut self, desired: Option<Token<F>>) -> Result<(), Error> {
    let tok = self.token_peek();
    if tok != desired {
      Err(Error::InvalidExpression)
    } else {
      self.token_next()?;
      Ok(())
    }
  }
}

fn eval_unary<F: Field>(lex: &mut Lexer<F>) -> Result<F, Error> {
  match lex.token_peek() {
    Some(Token::Number(num)) => {
      lex.token_next()?;
      Ok(num)
    },
    Some(Token::LParen) => {
      lex.token_next()?;
      let res = eval_expr(lex)?;
      lex.expect(Some(Token::RParen))?;
      Ok(res)
    },
    _ => Err(Error::InvalidExpression),
  }
}

fn eval_expr_prec<F: Field>(lex: &mut Lexer<F>, prec: usize) -> Result<F, Error> {
  if prec == 0 {
    return eval_unary(lex);
  }

  let mut result = eval_expr_prec(lex, prec-1)?;

  loop {
    let op = match lex.token_peek() {
      Some(Token::Oper(op)) if op.prec() == prec => {
        lex.token_next()?;
        op
      },
      _ => return Ok(result),
    };

    let rhs = eval_expr_prec(lex, prec-1)?;

    result = match op {
      Operator::Add => result + rhs,
      Operator::Sub => result - rhs,
      Operator::Mul => result * rhs,
      Operator::Div => (result / rhs).map_err(|_| Error::DivisionByZero)?,
    };
  }
}

fn eval_expr<F: Field>(lex: &mut Lexer<F>) -> Result<F, Error> {
  eval_expr_prec(lex, MAX_PREC)
}

pub fn evaluate<F: Field>(input: &str) -> Result<F, Error> {
  let mut lex = Lexer::new(input)?;
  let res = eval_expr(&mut lex)?;
  lex.expect(None)?;
  Ok(res)
}

//// A routine to pretty-print tables for a generic operation
fn print_operation_table<F: Field>(op: Operator) {
  let n = F::number_of_elements();

  let name = match op {
    Operator::Add => "+",
    Operator::Sub => "-",
    Operator::Mul => "*",
    Operator::Div => "/",
  };

  // Print the header
  print!(" {:>3}  | ", name);
  for x in 0..n {
    print!("{:>3}  ", x);
  }
  println!("");

  // Print the header separator
  println!("{}", "-".repeat(7 + 5*n));

  // Print each row
  for y in 0..n {
    let y_num: F = format!("{}", y).parse().unwrap();
    print!(" {:>3}  | ", y);
    // Print each column element
    for x in 0..n {
      let x_num: F = format!("{}", x).parse().unwrap();
      let result = match op {
        Operator::Add => Ok(y_num + x_num),
        Operator::Sub => Ok(y_num - x_num),
        Operator::Mul => Ok(y_num * x_num),
        Operator::Div => y_num / x_num,
      };
      let s = match result {
        Ok(num) => format!("{}", num),
        Err(_) => "-".to_string(),
      };
      print!("{:>3}  ", s);
    }
    println!("");
  }
}

fn read_line() -> String {
  let mut line = String::new();
  io::stdin().lock().read_line(&mut line).unwrap();
  line
}

pub fn interactive_calculator<F: Field>(name: &str) {
  let echo_line = match std::env::var("LEARN_YOU_GALOIS_FIELDS_FOR_GREAT_GOOD_CALC_ECHO") {
    Ok(val) => val != "0",
    Err(_) => false,
  };

  println!("{}", "=".repeat(80));
  println!("A Calculator for {}", name);
  println!("{}", "=".repeat(80));
  println!("");
  println!("Addition Table:");
  println!("");
  print_operation_table::<F>(Operator::Add);
  println!("");
  println!("Subtraction Table:");
  println!("");
  print_operation_table::<F>(Operator::Sub);
  println!("");
  println!("Multiplication Table:");
  println!("");
  print_operation_table::<F>(Operator::Mul);
  println!("");
  println!("Division Table:");
  println!("");
  print_operation_table::<F>(Operator::Div);
  println!("");
  println!("Enter any expression for evaluation (e.g. (1 + 2) * 4)");
  println!("");

  loop {
    print!("> ");
    io::stdout().flush().unwrap();

    let line = read_line();
    if line.len() == 0 {
      break;
    }
    if echo_line {
      print!("{}", line);
    }

    match evaluate::<F>(&line) {
      Ok(result) => println!("{}", result),
      Err(err) => println!("Error: {:?}", err),
    };
  }
}
