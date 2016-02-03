#include <vector>

namespace Check_pseudo_externals{
  template <std::size_t n>
  class Control_table{
  public:
    Control_table();
    Control_table(const Control_table&) = delete;
    Control_table(Control_table&&) = delete;
    Control_table& operator=(const Control_table&) = delete;
    Control_table& operator=(Control_table&&) = delete;
    ~Control_table() = default;

    template<std::size_t>
    bool not_init() const;
    template<std::size_t>
    void just_init();
  private:
    std::vector<bool> vb;
  }; 
}
template<std::size_t n>
Check_pseudo_externals::Control_table<n>::Control_table()
  : vb(n,false)
{
  static_assert(n>0, "Construct 'Control_table' with an positive integer.");
}
template<std::size_t n>
template<std::size_t i>
bool Check_pseudo_externals::Control_table<n>::not_init() const{
  static_assert(i < n, "Error in not_init().  "
		"The index must be between 0 and n-1");
  return !vb[i];
}
template<std::size_t n>
template<std::size_t i>
void Check_pseudo_externals::Control_table<n>::just_init(){
  static_assert(i < n, "Error in just_init().  "
		"The index must be between 0 and n-1");
  vb[i] = true;
}
