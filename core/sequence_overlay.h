#ifndef DELPHY_SEQUENCE_OVERLAY_H_
#define DELPHY_SEQUENCE_OVERLAY_H_

#include <absl/container/flat_hash_map.h>

#include "sequence.h"

namespace delphy {

// This whole structure is modelled on stl::vector<bool> from libstdc++

namespace details {

struct Sequence_overlay_impl {
  const Real_sequence* base_sequence_;
  absl::flat_hash_map<Site_index, Real_seq_letter> deltas_;
  
  auto get(Site_index l) const -> Real_seq_letter {
    if (auto it = deltas_.find(l); it != deltas_.end()) {
      return it->second;
    } else {
      return (*base_sequence_)[l];
    }
  }
  auto set(Site_index l, Real_seq_letter state) -> void {
    if ((*base_sequence_)[l] == state) {
      deltas_.erase(l);
    } else {
      deltas_[l] = state;
    }
  }
};

class Sequence_overlay_reference {
 public:
  Sequence_overlay_reference() : impl_{nullptr}, site_{0} {}
  Sequence_overlay_reference(Sequence_overlay_impl* impl, Site_index site) : impl_{impl}, site_{site} {}
  Sequence_overlay_reference(const Sequence_overlay_reference&) = default;

  operator Real_seq_letter() const {
    return impl_->get(site_); }
  auto operator=(Real_seq_letter state) -> Sequence_overlay_reference& {
    impl_->set(site_, state); return *this; }
  auto operator=(const Sequence_overlay_reference& x) noexcept -> Sequence_overlay_reference& {
    return *this = Real_seq_letter(x); }
  
  // See https://stackoverflow.com/a/63413121
  auto operator=(Real_seq_letter state) const -> void;  // <-- Not implemented, but signals AOK to std::indirectly_writable
  auto operator=(const Sequence_overlay_reference& x) const -> void;
  
  auto operator==(const Sequence_overlay_reference& x) const -> bool {
    return Real_seq_letter{*this} == Real_seq_letter{x}; };
  auto operator<=>(const Sequence_overlay_reference& that) const {
    return Real_seq_letter{*this} <=> Real_seq_letter{that}; }

  // Without the following stream writer, absl's CHECK_op macros get confused when dealing with these proxy references
  friend auto operator<<(std::ostream& os, const Sequence_overlay_reference& x) -> std::ostream& {
    return os << Real_seq_letter{x};
  }
  
 private:
  Sequence_overlay_impl* impl_;
  Site_index site_;
};

inline auto swap(Sequence_overlay_reference x, Sequence_overlay_reference y) noexcept -> void {
  auto tmp = Real_seq_letter{x};
  x = y;
  y = tmp;
}
inline auto swap(Sequence_overlay_reference x, Real_seq_letter& y) noexcept -> void {
  auto tmp = Real_seq_letter{x};
  x = y;
  y = tmp;
}
inline auto swap(Real_seq_letter& x, Sequence_overlay_reference y) noexcept -> void {
  auto tmp = Real_seq_letter{x};
  x = y;
  y = tmp;
}


class Sequence_overlay_iterator_base {
 public:
  using iterator_category = std::random_access_iterator_tag;
  using value_type = Real_seq_letter;
  using difference_type = Site_index;
  using pointer = void;

  Sequence_overlay_iterator_base(Sequence_overlay_impl* impl, Site_index site) : impl_{impl}, site_{site} {}
  
  friend auto operator==(
      const Sequence_overlay_iterator_base& x,
      const Sequence_overlay_iterator_base& y)
      -> bool {
    return x.impl_ == y.impl_ && x.site_ == y.site_;
  }
  friend auto operator<=>(
      const Sequence_overlay_iterator_base& x,
      const Sequence_overlay_iterator_base& y) noexcept {
    if (auto cmp = x.impl_ <=> y.impl_; cmp != 0) {
      return cmp;
    }
    return x.site_ <=> y.site_;
  }
  friend auto operator-(
      const Sequence_overlay_iterator_base& x,
      const Sequence_overlay_iterator_base& y)
      -> difference_type {
    return x.site_ - y.site_;
  }

 protected:
  Sequence_overlay_impl* impl_;
  Site_index site_;
};


class Sequence_overlay_iterator : public Sequence_overlay_iterator_base {
 public:
  using reference = Sequence_overlay_reference;
  using iterator = Sequence_overlay_iterator;

  Sequence_overlay_iterator()
      : Sequence_overlay_iterator_base{nullptr, 0} {}
  Sequence_overlay_iterator(Sequence_overlay_impl* impl, Site_index site)
      : Sequence_overlay_iterator_base{impl, site} {}

  auto operator*() const -> reference { return reference{impl_, site_}; }
  auto operator++() -> iterator& { ++site_; return *this; }
  auto operator++(int) -> iterator {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }
  auto operator--() -> iterator& { --site_; return *this; }
  auto operator--(int) -> iterator {
    auto tmp = *this;
    --(*this);
    return tmp;
  }

  auto operator+=(difference_type i) -> iterator& { site_ += i; return *this; }
  auto operator-=(difference_type i) -> iterator& { site_ -= i; return *this; }

  auto operator[](difference_type i) const -> reference { return *(*this + i); }

  friend auto operator+(const iterator& x, difference_type n) -> iterator {
    auto tmp = x;
    tmp += n;
    return tmp;
  }
  friend auto operator+(difference_type n, const iterator& x) -> iterator { return x + n; }
  friend auto operator-(const iterator& x, difference_type n) -> iterator {
    auto tmp = x;
    tmp -= n;
    return tmp;
  }
};

class Sequence_overlay_const_iterator : public Sequence_overlay_iterator_base {
 public:
  using reference = Real_seq_letter;
  using const_reference = Real_seq_letter;
  using const_iterator = Sequence_overlay_const_iterator;

  Sequence_overlay_const_iterator()
      : Sequence_overlay_iterator_base{nullptr, 0} {}
  Sequence_overlay_const_iterator(Sequence_overlay_impl* impl, Site_index site)
      : Sequence_overlay_iterator_base{impl, site} {}

  auto operator*() const -> const_reference { return impl_->get(site_); }
  auto operator++() -> const_iterator& { ++site_; return *this; }
  auto operator++(int) -> const_iterator {
    auto tmp = *this;
    ++(*this);
    return tmp;
  }
  auto operator--() -> const_iterator& { --site_; return *this; }
  auto operator--(int) -> const_iterator {
    auto tmp = *this;
    --(*this);
    return tmp;
  }

  auto operator+=(difference_type i) -> const_iterator& { site_ += i; return *this; }
  auto operator-=(difference_type i) -> const_iterator& { site_ -= i; return *this; }

  auto operator[](difference_type i) const -> const_reference { return *(*this + i); }

  friend auto operator+(const const_iterator& x, difference_type n) -> const_iterator {
    auto tmp = x;
    tmp += n;
    return tmp;
  }
  friend auto operator-(const const_iterator& x, difference_type n) -> const_iterator {
    auto tmp = x;
    tmp -= n;
    return tmp;
  }
  friend auto operator+(difference_type n, const const_iterator& x) -> const_iterator{ return x + n; }
};

}  // namespace details

// A non-owning view of a Real_sequence, with some modifications on top.
// Internally, keeps a reference to the base sequence and a map of sites to modified states
class Sequence_overlay {
 public:
  using value_type = Real_sequence::value_type;
  using size_type = Site_index;
  using difference_type = Site_index;
  using reference = details::Sequence_overlay_reference;
  using const_reference = Real_seq_letter;
  using pointer = void;
  using const_pointer = void;
  using iterator = details::Sequence_overlay_iterator;
  using const_iterator = details::Sequence_overlay_const_iterator;
  
  explicit Sequence_overlay(const Real_sequence& base_sequence,
                            absl::flat_hash_map<Site_index, Real_seq_letter> deltas = {})
      : impl_{std::addressof(base_sequence), std::move(deltas)} {}

  auto base_sequence() const -> const Real_sequence& { return *impl_.base_sequence_; }
  auto deltas() const -> const absl::flat_hash_map<Site_index, Real_seq_letter> { return impl_.deltas_; }

  auto begin() noexcept -> iterator {
    return iterator{&impl_, 0}; }
  auto begin() const noexcept -> const_iterator {
    return const_iterator{const_cast<details::Sequence_overlay_impl*>(&impl_), 0}; }
  auto end() noexcept -> iterator {
    return iterator{&impl_, size()}; }
  auto end() const noexcept -> const_iterator {
    return const_iterator{const_cast<details::Sequence_overlay_impl*>(&impl_), size()}; }
  
  auto size() const -> size_type { return static_cast<Site_index>(impl_.base_sequence_->size()); }
  auto empty() const -> bool { return impl_.base_sequence_->empty(); }

  auto operator[](difference_type l) -> reference { return {&impl_, l}; }
  auto operator[](difference_type l) const -> const_reference { return impl_.get(l); }

  auto materialize() const -> Real_sequence {
    auto result = Real_sequence{*impl_.base_sequence_};  // copy
    for (auto [l, s] : impl_.deltas_) {
      result[l] = s;
    }
    return result;
  }

  auto rebase_to(const Real_sequence& new_base_sequence) -> Sequence_overlay& {
    const auto& old_base_sequence = *impl_.base_sequence_;
    CHECK_EQ(std::ssize(old_base_sequence), std::ssize(new_base_sequence));
    
    impl_.base_sequence_ = std::addressof(new_base_sequence);

    for (auto l = Site_index{0}; l != std::ssize(new_base_sequence); ++l) {
      if (old_base_sequence[l] != new_base_sequence[l]) {
        if (auto it = impl_.deltas_.find(l); it != impl_.deltas_.end()) {
          // X -> Z should become Y -> Z
          auto new_from = new_base_sequence[l];
          auto to = it->second;
          if (new_from == to) {
            // if Y == Z, then the delta disappears
            impl_.deltas_.erase(it);
          }
        } else {
          // X -> X should become Y -> X
          impl_.deltas_.insert({l, old_base_sequence[l]});
        }
      }
    }

    return *this;
  }

 private:
  details::Sequence_overlay_impl impl_;
};

}  // namespace delphy

#endif // DELPHY_SEQUENCE_OVERLAY_H_
