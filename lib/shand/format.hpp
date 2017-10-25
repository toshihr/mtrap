#ifndef SHAND_FORMAT_INCLUDE
#define SHAND_FORMAT_INCLUDE

/*
 http://d.hatena.ne.jp/faith_and_brave/20080417/1208430194

using namespace std;
using namespace shand;
int main()
{
    // printf風
    string str1 = format("Hello %05d %s %f") % 3 % "abc" % 3.14;

    // プレースホルダー(Boost.Format風)
    string str2 = format("Hello %1% %2% %3%") % 3 % "abc" % 3.14;

    // プレースホルダー(C#風)
    string str3 = format("Hello {0} {1} {2}") % 3 % "abc" % 3.14;

    // wstring版
    wstring wstr = wformat(L"Hello %05d %s %f") % 3 % L"abc" % 3.14;

    return 0;
}
*/

#include <iomanip>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm> // find

namespace shand {

namespace format_detail {
    template <class CharT, class T>
    inline std::basic_string<CharT> to_string(const T& x)
    {
        std::basic_stringstream<CharT> interpreter;
        interpreter << x;
        return interpreter.str();
    }

    template <class StringT>
    inline StringT replace(StringT source, const StringT& oldstr, const StringT& newstr)
    {
        size_t i;
        while ((i = source.find(oldstr)) != StringT::npos)
            source.replace(i, oldstr.length(), newstr);
        return source;
    }

} // namespace format_detail

// フォーマットの解析
template <class CharT, class Traits=std::char_traits<CharT> >
class printf_formatter {
    typedef CharT                                       char_type;
    typedef std::basic_string<char_type, Traits>        string_type;
    typedef std::basic_stringstream<char_type, Traits>  stream_type;
    typedef typename string_type::const_iterator        string_iterator;

    const string_type&              fmt_;   // フォーマット
    const std::vector<string_type>& items_; // 要素
    std::vector<char_type>          types_; // フォーマットのコード

public:
    typedef string_type result_type;

    printf_formatter(const string_type& fmt, const std::vector<string_type>& items)
        : fmt_(fmt), items_(items)
    {
        types_.reserve(14);
        types_.push_back('d'); types_.push_back('i'); types_.push_back('o');
        types_.push_back('u'); types_.push_back('x'); types_.push_back('X');
        types_.push_back('e'); types_.push_back('E'); types_.push_back('f');
        types_.push_back('F'); types_.push_back('g'); types_.push_back('G');
        types_.push_back('c'); types_.push_back('s');
    }

    result_type operator()() const
    {
        string_type result;
        const char_type rep      = '%';
        string_iterator pos      = fmt_.begin();
        string_iterator fmt_last = fmt_.end();

        // %が見つかるまで無視
        while (pos != fmt_last && *pos != rep)
            result += *pos++;

        typedef typename std::vector<string_type>::const_iterator item_iterator;
        for (item_iterator it = items_.begin(), last = items_.end(); it != last; ++it) {
            if (*pos == rep && pos + 1 != fmt_last && *(pos + 1) != rep) {
                string_type flag;
                while (++pos != fmt_last && is_flag(*pos))
                    flag += *pos;

                if (pos == fmt_last)
                    return result;

                // フォーマット(%??)を置換
                result += put(*it, *pos++, flag);
            }

            // %が見つかるまで無視
            while (pos != fmt_last && *pos != rep)
                result += *pos++;
        }

        return result;
    }

    bool is_flag(char_type c) const
    {
        return std::find(types_.begin(), types_.end(), c) == types_.end();
    }

    template <class T>
    string_type put(const T& arg, char_type type, const string_type& flag) const
    {
        stream_type ss;

        flags_field(ss, flag);
        type_char_field(ss, type);
        ss << arg;

        return ss.str();
    }

    void flags_field(stream_type& ss, const string_type& flag) const
    {
        const char_type space = ' ';
        bool is_left = false;

        for (string_iterator pos = flag.begin(), last = flag.end(); pos != last; ++pos) {
            switch (*pos) {
            // 左寄せ
            case '-':
                is_left = true;
                ss << std::setfill(space);
                ss << std::left;
                break;

            // 8/16進数に接頭辞を出力
            case '#':
                ss << std::showbase;
                break;

            // 正数に+を出力
            case '+':
                ss << std::showpos;
                break;

            // 余白
            case ' ':
                ss << std::setfill(*pos);
                break;

            // 0埋め
            case '0':
                if (!is_left)
                    ss << std::setfill(*pos);
                break;

            // 小数
            case '.':
                set_precision(ss, ++pos, last);
                if (pos == last)
                    return;
                break;

            default:
                if (isdigit(*pos))
                    set_width(ss, pos, last);

                if (pos == last)
                    return;
                break;
            }
        }
    }

    int to_value(const string_type& str) const
    {
        int value = 0;
        stream_type interpreter;
        interpreter << str;
        interpreter >> value;
        return value;
    }

    void set_width(stream_type& ss, string_iterator& cur, string_iterator last) const
    {
        string_type digit;
        for (; cur != last && isdigit(*cur); ++cur)
            digit += *cur;

        if (!digit.empty())
            ss << std::setw(to_value(digit));

        if (cur != last && *cur == '.')
            set_precision(ss, ++cur, last);
    }

    void set_precision(stream_type& ss, string_iterator& cur, string_iterator last) const
    {
        string_type digit;
        for (; cur != last && isdigit(*cur); ++cur)
            digit += *cur;

        if (!digit.empty())
            ss << std::setprecision(to_value(digit));
    }

    void type_char_field(stream_type& ss, char_type type) const
    {
        if (isupper(type))
            ss << std::uppercase;

        switch (tolower(type)) {
        // 8進数
        case 'o':
            ss << std::oct;
            break;

        // 16進数
        case 'x':
            ss << std::hex;
            break;

        // 浮動小数点(指数表記)
        case 'e':
            ss << std::scientific;
            break;

        // 浮動小数点(固定小数点表記)
        case 'f':
            ss << std::setiosflags(std::ios::fixed);
            break;

        default:
            break;
        }
    }
};

// プレースホルダーの解析
template <class CharT, class Traits=std::char_traits<CharT> >
class placeholder_formatter {
    typedef std::basic_string<CharT, Traits>        string_type;
    typedef std::basic_stringstream<CharT, Traits>  stream_type;

    const string_type&              fmt_;
    const std::vector<string_type>& items_;
public:
    typedef string_type result_type;

    placeholder_formatter(const string_type& fmt, const std::vector<string_type>& items)
        : fmt_(fmt), items_(items) {}

    result_type operator()() const
    {
        string_type result = fmt_;
        for (size_t i = 0; i < items_.size(); ++i) {
            // Boost.Format風
            {
                stream_type ss;
                ss << "%";
                ss << i + 1;
                ss << "%";
                result = format_detail::replace<string_type>(result, ss.str(), items_[i]);
            }

            // C#風
            {
                stream_type ss;
                ss << "{";
                ss << i;
                ss << "}";
                result = format_detail::replace<string_type>(result, ss.str(), items_[i]);
            }
        }

        return result;
    }
};

// basic_format
template <class CharT, class Traits=std::char_traits<CharT> >
class basic_format {
public:
    typedef std::basic_string<CharT, Traits> string_type;

    basic_format(const string_type& fmt)
        : fmt_(fmt) {}

    template <class T>
    basic_format& operator%(const T& arg)
    {
        items_.push_back(format_detail::to_string<CharT>(arg));
        return *this;
    }

    operator string_type() const { return analysis(); }
    string_type str() const { return analysis(); }

private:
    string_type analysis() const
    {
        string_type result;
        result = placeholder_formatter<CharT>(fmt_, items_)();  // プレースホルダーの解析
        if (result == fmt_)
            result = printf_formatter<CharT>(fmt_, items_)();   // フォーマットの解析
        return result;
    }

    string_type fmt_;
    std::vector<string_type> items_;
};

typedef basic_format<char>      format;
typedef basic_format<wchar_t>   wformat;

} // namespace shand


#endif // SHAND_FORMAT_INCLUDE

