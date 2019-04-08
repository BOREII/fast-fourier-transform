#include "substring_matching.h"
#include "test_runner.h"
#include "profile.h"

namespace SubstringMatching {

void TestSubstrings() {
    std::string s = "aaaa";
    ASSERT_EQUAL(FindSubstrings(s, "a"), FindSubstringsFFT(s, "a"));
    ASSERT_EQUAL(FindSubstringsFFT("ababab", "ababab"), FindSubstrings("ababab", "ababab"))
}

void TestSubstringsHeyJude() {
    std::string hey_jude_lyrics = R"(Hey, Jude, don't make it bad
Take a sad song and make it better
Remember to let her into your heart
Then you can start to make it better

Hey, Jude, don't be afraid
You were made to go out and get her
The minute you let her under your skin
Then you begin to make it better

And anytime you feel the pain,
Hey, Jude, refrain
Don't carry the world upon your shoulders
For well you know that it's a fool
Who plays it cool
By making his world a little colder

Nah, nah nah, nah nah, nah nah, nah nah

Hey, Jude, don't let me down
You have found her, now go and get her
Remember to let her into your heart
Then you can start to make it better

So let it out and let it in,
Hey, Jude, begin
You're waiting for someone to perform with
And don't you know that it's just you,
Hey, Jude, you'll do
The movement you need is on your shoulder

Nah, nah nah, nah nah, nah nah, nah nah yeah

Hey, Jude, don't make it bad
Take a sad song and make it better
Remember to let her under your skin
Then you'll begin to make it better, better, better, better, better... oh!

Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (Jude)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (yeah, yeah, yeah)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (don't make it bad, Jude)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (take a sad song and make it better)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (oh, Jude)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (Jude, hey, Jude, whoa)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (ooh)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude)";

    std::string word_Hey = "Hey", word_Jude = "Jude", word_nah = "nah";

    std::vector<size_t> Hey_positions, Jude_positions, nah_positions,
                        Hey_positions_fft, Jude_positions_fft, nah_positions_fft;

    {
        LOG_DURATION("Substrings:n^2 algorithm");
        Hey_positions = FindSubstrings(hey_jude_lyrics, word_Hey);
        Jude_positions = FindSubstrings(hey_jude_lyrics, word_Jude);
        nah_positions = FindSubstrings(hey_jude_lyrics, word_nah);

    }
    {
        LOG_DURATION("Substrings:nlogn algorithm");
        Hey_positions_fft = FindSubstringsFFT(hey_jude_lyrics, word_Hey);
        Jude_positions_fft = FindSubstringsFFT(hey_jude_lyrics, word_Jude);
        nah_positions_fft = FindSubstringsFFT(hey_jude_lyrics, word_nah);
    }

    ASSERT_EQUAL(Hey_positions, Hey_positions_fft);
    ASSERT_EQUAL(Jude_positions, Jude_positions_fft);
    ASSERT_EQUAL(nah_positions, nah_positions_fft);
}

void TestMatches() {
    ASSERT_EQUAL(FindMatches("aaaa", "a?"), FindMatchesFFT("aaaa", "a?"));
    ASSERT_EQUAL(FindMatches("ababab", "a??b?b"), FindMatchesFFT("ababab", "a??b?b"))
}

void TestMatchesHeyJude() {
    std::string hey_jude_lyrics = R"(Hey, Jude, don't make it bad
Take a sad song and make it better
Remember to let her into your heart
Then you can start to make it better

Hey, Jude, don't be afraid
You were made to go out and get her
The minute you let her under your skin
Then you begin to make it better

And anytime you feel the pain,
Hey, Jude, refrain
Don't carry the world upon your shoulders
For well you know that it's a fool
Who plays it cool
By making his world a little colder

Nah, nah nah, nah nah, nah nah, nah nah

Hey, Jude, don't let me down
You have found her, now go and get her
Remember to let her into your heart
Then you can start to make it better

So let it out and let it in,
Hey, Jude, begin
You're waiting for someone to perform with
And don't you know that it's just you,
Hey, Jude, you'll do
The movement you need is on your shoulder

Nah, nah nah, nah nah, nah nah, nah nah yeah

Hey, Jude, don't make it bad
Take a sad song and make it better
Remember to let her under your skin
Then you'll begin to make it better, better, better, better, better... oh!

Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (Jude)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (yeah, yeah, yeah)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (don't make it bad, Jude)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (take a sad song and make it better)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (oh, Jude)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (Jude, hey, Jude, whoa)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude (ooh)
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude
Nah, nah nah, nah nah, nah, nah, nah nah,
Hey, Jude)";

    std::string word_Hey = "?e?", word_Jude = "J?de", word_nah = "?ah";

    std::vector<size_t> Hey_positions, Jude_positions, nah_positions,
        Hey_positions_fft, Jude_positions_fft, nah_positions_fft;

    {
        LOG_DURATION("Matches:n^2 algorithm");
        Hey_positions = FindMatches(hey_jude_lyrics, word_Hey);
        Jude_positions = FindMatches(hey_jude_lyrics, word_Jude);
        nah_positions = FindMatches(hey_jude_lyrics, word_nah);

    }
    {
        LOG_DURATION("Matches:nlogn algorithm");
        Hey_positions_fft = FindMatchesFFT(hey_jude_lyrics, word_Hey);
        Jude_positions_fft = FindMatchesFFT(hey_jude_lyrics, word_Jude);
        nah_positions_fft = FindMatchesFFT(hey_jude_lyrics, word_nah);
    }

    ASSERT_EQUAL(Hey_positions, Hey_positions_fft);
    ASSERT_EQUAL(Jude_positions, Jude_positions_fft);
    ASSERT_EQUAL(nah_positions, nah_positions_fft);
}
}