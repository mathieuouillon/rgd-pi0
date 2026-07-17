#pragma once

/// \file Sha256.hpp
/// \brief SHA-256 (FIPS 180-4) over a byte range or a file. Header-only.
///
/// WHY THIS EXISTS AS A HEADER RATHER THAN IN THE ONE PROGRAM THAT NEEDED IT
/// FIRST: this hash is a PROVENANCE FINGERPRINT, and its whole value is that
/// the same bytes give the same string in every stage. stageA_skim stamps
/// `config.sha256` into its output; make_grid stamps `config.sha256` into the
/// grid JSON and COMPARES it against the one in the slim files it reads. Those
/// two numbers are only comparable if they come from one implementation. Two
/// copies of a hash function are two copies of a constant: they agree until the
/// day one is touched, and then every downstream cross-check silently reports a
/// mismatch that is an artefact of the code rather than of the config.
///
/// Hand-rolled rather than pulled from OpenSSL: the project has no crypto
/// dependency and this is a fingerprint, not a security boundary. Nothing here
/// resists a chosen-prefix attack and nothing needs to -- the threat model is
/// "did somebody edit cuts.json between the skim and the binning", not an
/// adversary.
///
/// Verified against `shasum -a 256` on config/cuts.json.
///
/// Header-only (`inline`) so it needs no meson target of its own: every
/// consumer already carries `src/` as an include directory via one of the
/// pi0_*_dep dependencies. If it ever grows past a page or acquires state worth
/// testing on its own, promote it to src/util/Sha256.cpp and a pi0_util
/// library; there is nothing here that makes that hard.

#include <array>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace pi0::util {

/// Streaming SHA-256. Feed it with update(), read it once with final_hex().
///
/// final_hex() applies the padding and is therefore NOT idempotent: calling it
/// twice on one object gives a second answer that is not the digest of
/// anything. One object, one digest.
class Sha256 {
   public:
    void update(const std::uint8_t* data, std::size_t len) {
        for (std::size_t i = 0; i < len; ++i) {
            m_buf[m_buflen++] = data[i];
            if (m_buflen == 64) {
                transform();
                m_bitlen += 512;
                m_buflen = 0;
            }
        }
    }

    [[nodiscard]] std::string final_hex() {
        std::size_t i = m_buflen;
        m_buf[i++] = 0x80;  // the mandatory trailing 1 bit
        if (m_buflen >= 56) {
            while (i < 64) m_buf[i++] = 0x00;
            transform();
            i = 0;
        }
        while (i < 56) m_buf[i++] = 0x00;

        m_bitlen += static_cast<std::uint64_t>(m_buflen) * 8;
        for (int b = 0; b < 8; ++b) m_buf[63 - static_cast<std::size_t>(b)] = static_cast<std::uint8_t>(m_bitlen >> (8 * b));
        transform();

        std::ostringstream os;
        os << std::hex << std::setfill('0');
        for (std::uint32_t h : m_state) os << std::setw(8) << h;
        return os.str();
    }

   private:
    static std::uint32_t rotr(std::uint32_t x, std::uint32_t n) { return (x >> n) | (x << (32 - n)); }

    void transform() {
        static constexpr std::array<std::uint32_t, 64> kK = {
            0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
            0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
            0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
            0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
            0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
            0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
            0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
            0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
        };

        std::array<std::uint32_t, 64> w{};
        for (std::size_t i = 0; i < 16; ++i) {
            w[i] = (static_cast<std::uint32_t>(m_buf[i * 4]) << 24) | (static_cast<std::uint32_t>(m_buf[i * 4 + 1]) << 16) |
                   (static_cast<std::uint32_t>(m_buf[i * 4 + 2]) << 8) | static_cast<std::uint32_t>(m_buf[i * 4 + 3]);
        }
        for (std::size_t i = 16; i < 64; ++i) {
            const std::uint32_t s0 = rotr(w[i - 15], 7) ^ rotr(w[i - 15], 18) ^ (w[i - 15] >> 3);
            const std::uint32_t s1 = rotr(w[i - 2], 17) ^ rotr(w[i - 2], 19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16] + s0 + w[i - 7] + s1;
        }

        std::uint32_t a = m_state[0], b = m_state[1], c = m_state[2], d = m_state[3];
        std::uint32_t e = m_state[4], f = m_state[5], g = m_state[6], h = m_state[7];

        for (std::size_t i = 0; i < 64; ++i) {
            const std::uint32_t s1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            const std::uint32_t ch = (e & f) ^ (~e & g);
            const std::uint32_t t1 = h + s1 + ch + kK[i] + w[i];
            const std::uint32_t s0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            const std::uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
            const std::uint32_t t2 = s0 + maj;
            h = g; g = f; f = e; e = d + t1;
            d = c; c = b; b = a; a = t1 + t2;
        }

        m_state[0] += a; m_state[1] += b; m_state[2] += c; m_state[3] += d;
        m_state[4] += e; m_state[5] += f; m_state[6] += g; m_state[7] += h;
    }

    std::array<std::uint32_t, 8> m_state = {
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19,
    };
    std::array<std::uint8_t, 64> m_buf{};
    std::uint64_t m_bitlen = 0;
    std::size_t m_buflen = 0;
};

/// SHA-256 of a file's bytes, lowercase hex.
/// \throws std::runtime_error if the file cannot be read -- a provenance record
///         with a blank hash is worse than no file at all.
[[nodiscard]] inline std::string sha256_file(const std::string& path) {
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::runtime_error("cannot open for hashing: " + path);

    Sha256 h;
    std::array<char, 64 * 1024> buf{};
    while (in.read(buf.data(), static_cast<std::streamsize>(buf.size())) || in.gcount() > 0) {
        h.update(reinterpret_cast<const std::uint8_t*>(buf.data()), static_cast<std::size_t>(in.gcount()));
        if (!in) break;
    }
    return h.final_hex();
}

}  // namespace pi0::util
