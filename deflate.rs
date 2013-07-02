extern mod std;
extern mod extra;

use std::io;
use std::io::{Reader,ReaderUtil};
use std::uint;
use extra::io_util::BufReader;

static max_bits:uint = 15;
static max_lit_codes:uint = 286;
static max_dist_codes:uint = 30;
static max_codes:uint = max_lit_codes + max_dist_codes;
static fixed_lit_codes:uint = 288;

struct HuffmanTable {
    counts: ~[uint],
    symbols: ~[uint],
}

fn huffman_new(sym_lens: &[uint]) -> Option<~HuffmanTable> {
    let len = sym_lens.len();
    let mut counts = ~[];
    counts.grow(len, &0);
    for sym_lens.iter().advance |len| {
        counts[*len] += 1;
    }

    let mut left = 1;
    for uint::range(1, max_bits+1) |len| {
        left <<= 1;
        left -= counts[len];
        if(left < 0) {
            return None;
        }
    }

    let mut offs = ~[];
    offs.grow(max_bits+1, &0);
    for uint::range(1, max_bits) |len| {
        offs[len+1] = offs[len] + counts[len];
    }

    let mut symbols = ~[];
    symbols.grow(len, &0);
    for uint::range(0, len) |i| {
        if(sym_lens[i] != 0) {
            let offs_index = sym_lens[i];
            symbols[offs[offs_index]] = i;
            offs[offs_index] += 1;
        }
    }

    Some(~HuffmanTable{ counts: counts, symbols: symbols })
}

struct DeflateStream {
    bitPos: uint,
    bitBuf: uint,

    reader: @Reader,
    writeBuf: ~[u8],
}

fn new(reader: @Reader) -> ~DeflateStream {
    return ~DeflateStream {
        bitPos: 0,
        bitBuf: 0,
        reader: reader,
        writeBuf: ~[],
    }
}

static order: &'static[u8] = &'static[
            16, 17, 18, 0, 8, 7, 9, 6,
            10, 5, 11, 4, 12, 3, 13, 2,
            14, 1, 15];
static lookup_lens: &'static[uint] = &'static[
            3, 4, 5, 6, 7, 8, 9, 10,
            11, 13, 15, 17, 19, 23, 27, 31,
            35, 43, 51, 59, 67, 83, 99, 115,
            131, 163, 195, 227, 258];
static lookup_len_bits: &'static[uint] = &'static[
            0, 0, 0, 0, 0, 0, 0, 0,
            1, 1, 1, 1, 2, 2, 2, 2,
            3, 3, 3, 3, 4, 4, 4, 4,
            5, 5, 5, 5, 0];
static lookup_dists: &'static[uint] = &'static[
            1, 2, 3, 4, 5, 7, 9, 13,
            17, 25, 33, 49, 65, 97, 129, 193,
            257, 385, 513, 769, 1025, 1537, 2049, 3073,
            4097, 6145, 8193, 12289, 16385, 24577];
static lookup_dist_bits: &'static[uint] = &'static[
            0, 0, 0, 0, 1, 1, 2, 2,
            3, 3, 4, 4, 5, 5, 6, 6,
            7, 7, 8, 8, 9, 9, 10, 10,
            11, 11, 12, 12, 13, 13];

impl DeflateStream {
    fn bits(&mut self, len: uint) -> Option<uint> {
        while(self.bitPos < len) {
            let next = self.reader.read_byte();
            if(next < 0) {
                return None;
            }
            self.bitBuf |= (next as uint) << self.bitPos;
            self.bitPos += 8;
        }
        let val = self.bitBuf & ((1 << len) - 1);
        self.bitBuf = self.bitBuf >> len;
        self.bitPos -= len;
        Some(val)
    }

    fn parse_symbol(&mut self, table: &HuffmanTable) -> Option<uint> {
        let mut code = 0;
        let mut first = 0;
        let mut index = 0;
        for uint::range(1, max_bits+1) |len| {
            match self.bits(1) {
                Some(val) => {
                    code = (code << 1) | val;
                    let count = table.counts[len];
                    if(code < first + count) {
                        return Some(table.symbols[index + (code - first)])
                    }
                    index += count;
                    first = (first + count) << 1;
                },
                None => { return None; }
            }
        }
        None
    }

    fn decode(&mut self, len_table: &HuffmanTable, dist_table: &HuffmanTable) -> bool {
        loop {
            match(self.parse_symbol(len_table)) {
                Some(sym) if sym < 256 => {
                    self.writeBuf.push(sym as u8);
                },
                Some(sym) if sym == 256 => {
                    return true;
                },
                Some(sym) if sym < 257 + 29 => {
                    sym -= 257;
                    match(self.bits(lookup_len_bits[sym]), self.parse_symbol(dist_table)) {
                        (Some(bits), Some(dist_sym)) => {
                            match(self.bits(lookup_dist_bits[dist_sym])) {
                                Some(dist_bits) => {
                                    let len = lookup_lens[sym] + bits;
                                    let dist = lookup_dists[dist_sym] + dist_bits;
                                    if(self.writeBuf.len() < dist) {
                                        return false;
                                    }
                                    let base_idx = self.writeBuf.len() - dist;
                                    for uint::range(0, len) |i| {
                                        self.writeBuf.push(self.writeBuf[base_idx + i]);
                                    }
                                },
                                None => { return false; },
                            };
                        },
                        _ => { return false; }
                    };
                },
                _ => { return false; },
            };
        }
    }

    fn fixed(&mut self) -> bool {
        let mut lengths = ~[];
        for uint::range(0, 144) |_| {
            lengths.push(8);
        }
        for uint::range(144, 256) |_| {
            lengths.push(9);
        }
        for uint::range(256, 280) |_| {
            lengths.push(7);
        }
        for uint::range(280, fixed_lit_codes) |_| {
            lengths.push(8);
        }

        let mut dists = ~[];
        dists.grow(30, &5);

        match (huffman_new(lengths), huffman_new(dists)) {
            (Some(len_table), Some(dist_table)) => {
                self.decode(len_table, dist_table)
            }
            _ => false
        }
    }

    fn dyn_table(&mut self, code_len_table: &HuffmanTable, len: uint)
            -> Option<~HuffmanTable> {
        // parse literal table
        let mut lengths = ~[];
        lengths.grow(len, &0);

        let mut i = 0;
        while (i < len) {
            match(self.parse_symbol(code_len_table)) {
                Some(sym) if sym < 16 => {
                    lengths[i] = sym;
                    i += 1;
                },
                Some(sym) => {
                    let (sym, base_len, bit_len) = match(sym) {
                        16 => (lengths[i-1], 3, 2),
                        17 => (0, 3, 3),
                        18 => (0, 11, 7),
                        _ => { return None; }
                    };

                    match(self.bits(bit_len)) {
                        Some(len) => {
                            for uint::range(0, base_len + len) |_| {
                                lengths[i] = sym;
                                i += 1;
                            }
                        },
                        None => { return None; }
                    }
                },
                None => { return None; }
            }
        };
        huffman_new(lengths)
    }

    fn dynamic(&mut self) -> bool {
        match (self.bits(5), self.bits(5), self.bits(4)) {
            (Some(lenbits), Some(distbits), Some(codebits)) => {
                let nlen = lenbits + 257;
                let ndist = distbits + 1;
                let ncode = codebits + 4;
                if(nlen > max_lit_codes || ndist > max_dist_codes) {
                    return false;
                }

                let mut lengths = ~[];
                lengths.grow(19, &0);
                for uint::range(0, ncode) |i| {
                    match(self.bits(3)) {
                        Some(len) => {
                            lengths[order[i]] = len;
                        },
                        None => { return false; }
                    }
                }

                match(huffman_new(lengths)) {
                    Some(code_len_table) => {
                        match(self.dyn_table(code_len_table, nlen), 
                                self.dyn_table(code_len_table, ndist)) {
                            (Some(len_table), Some(dist_table)) => {
                                self.decode(len_table, dist_table)
                            },
                            _ => false
                        }
                    },
                    None => false
                }
            },
            _ => {
                false
            },
        }
    }

    fn inflate(&mut self) -> bool {
        loop {
            let (success, done) = match (self.bits(1), self.bits(2)) {
                (Some(last), Some(comp_type)) => {
                    match(comp_type) {
                        0 => {
                            fail!("Not implemented: store")
                        },
                        1 => {
                            (self.fixed(), last != 0)
                        },
                        2 => {
                            (self.dynamic(), last != 0)
                        },
                        _ => {
                            (false, false)
                        },
                    }
                },
                _ => {
                    (false, false)
                },
            };

            if(!success) {
                return false;
            }
            if(done) {
                return true;
            }
        }
    }
}

fn main() {
    let bytes = io::stdin().read_whole_stream();
    let bufReader = @BufReader::new(bytes);
    let mut decoder = new(bufReader as @Reader);
    if(decoder.inflate()) {
        print(fmt!("Success: %u bytes", decoder.writeBuf.len()));
    } else {
        print("Failed");
    }
}
