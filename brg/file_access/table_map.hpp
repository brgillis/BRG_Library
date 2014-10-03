/**********************************************************************\
 @file table_map.hpp
 ------------------

 TODO <Insert file description here>

 **********************************************************************

 Copyright (C) 2014  Bryan R. Gillis

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

\**********************************************************************/


// body file: table_map.cpp


#ifndef _BRG_TABLE_MAP_HPP_INCLUDED_
#define _BRG_TABLE_MAP_HPP_INCLUDED_

#include <functional>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <unordered_map>

namespace brgastro {

template<typename Key>
class compare_insertion_order
{
private:
	mutable const std::unordered_map<Key,size_t> *_key_map_;
public:

	compare_insertion_order()
	: _key_map_(NULL)
	{
	}
	compare_insertion_order(const std::unordered_map<Key,size_t> *new_key_map)
	: _key_map_(new_key_map)
	{
	}
	compare_insertion_order(const std::unordered_map<Key,size_t> &new_key_map)
	: _key_map_(&new_key_map)
	{
	}

	bool operator()(const Key & key1, const Key & key2) const
	{
		return _key_map_->at(key1) < _key_map_->at(key2);
	}

	void set_key_map(const std::unordered_map<Key,size_t> *new_key_map) const
	{
		_key_map_ = new_key_map;
	}
	void set_key_map(const std::unordered_map<Key,size_t> &new_key_map) const
	{
		_key_map_ = std::addressof(new_key_map);
	}
};

template <typename T>
class table_map_t
{
private:
	typedef std::map<std::string,std::vector<T>,std::function<bool(const std::string&,const std::string&)>> base_type;
	base_type base;

	mutable std::unordered_map<typename base_type::key_type,size_t> _key_map_;
	mutable size_t _counter_;

	mutable compare_insertion_order<typename base_type::key_type> _comparer_;

	void _update_key_map() const
	{
		_comparer_.set_key_map(_key_map_);
	}
	void _confirm_key(const typename base_type::key_type & k) const
	{
		if(_key_map_.count(k)==0) _key_map_[k]=_counter_++;
	}

public:

	// Overloads of the constructor to set up the key comparison operator properly
#if(1)
	table_map_t()
	: _counter_(0)
	{
		auto func = [&] (const typename base_type::key_type & key1, const typename base_type::key_type & key2)
				{return _comparer_(key1,key2);};
		base = base_type(func,typename base_type::allocator_type());
		_update_key_map();
	}
	explicit table_map_t(const typename base_type::allocator_type& alloc)
	: _counter_(0)
	{
		auto func = [&] (const typename base_type::key_type & key1, const typename base_type::key_type & key2)
				{return _comparer_(key1,key2);};
		base = base_type(func,alloc);
		_update_key_map();
	}
	template <class InputIterator>
	table_map_t (InputIterator first, InputIterator last,
		const typename base_type::allocator_type& alloc = typename base_type::allocator_type())
	: _counter_(0)
	{
		auto func = [&] (const typename base_type::key_type & key1, const typename base_type::key_type & key2)
				{return _comparer_(key1,key2);};
		base = base_type(first,last,func,alloc);
		_update_key_map();
    }
#endif

	// Overloads of copy, assignment, and move
#if(1)
	table_map_t(const table_map_t & other)
	: base(other.base),
	  _key_map_(other._key_map_),
	  _counter_(other._counter_),
	  _comparer_(other._comparer_)
	{
		_update_key_map();
	}
	table_map_t(table_map_t&& other)
	: table_map_t()
	{
		swap(other);
		_update_key_map();
	}
	table_map_t & operator=(const table_map_t & other)
	{
		swap(table_map_t(other));
		_update_key_map();
		return *this;
	}
	table_map_t & operator=(table_map_t&& other)
	{
		swap(other);
		_update_key_map();
		return *this;
	}
	virtual ~table_map_t() { }
#endif

#if(1) // overloads of map functions to update key map first

	typename base_type::mapped_type& at (const typename base_type::key_type& k)
	{
		return base.at(k);
	}
	const typename base_type::mapped_type& at (const typename base_type::key_type& k) const
	{
		return base.at(k);
	}

	typename base_type::iterator find (const typename base_type::key_type& k)
    {
    	if(count(k)==0) return base.end();
		return base.find(k);
    }
	typename base_type::const_iterator find (const typename base_type::key_type& k) const
    {
    	if(count(k)==0) return base.end();
		return base.find(k);
    }

    typename base_type::size_type count (const typename base_type::key_type& k) const
    {
		if (_key_map_.count(k)==0) return 0;
		return base.count(k);
    }

    typename base_type::iterator lower_bound (const typename base_type::key_type& k)
    {
    	if(count(k)==0) return base.end();
		return base.lower_bound(k);
    }
    typename base_type::const_iterator lower_bound (const typename base_type::key_type& k) const
    {
    	if(count(k)==0) return base.end();
		return base.lower_bound(k);
    }

    typename base_type::iterator upper_bound (const typename base_type::key_type& k)
    {
    	if(count(k)==0) return base.end();
		return base.upper_bound(k);
    }
    typename base_type::const_iterator upper_bound (const typename base_type::key_type& k) const
    {
    	if(count(k)==0) return base.end();
		return base.upper_bound(k);
    }

    std::pair<typename base_type::const_iterator,typename base_type::const_iterator> equal_range
    	(const typename base_type::key_type& k) const
	{
    	_confirm_key(k);
		return base.equal_range(k);
	}
    std::pair<typename base_type::iterator,typename base_type::iterator> equal_range (const typename base_type::key_type& k)
	{
    	_confirm_key(k);
		return base.equal_range(k);
	}

#endif

    // More complicated overloads
#if(1)

    typename base_type::mapped_type & operator[](const typename base_type::key_type & key)
	{
    	_confirm_key(key);
		return base.operator[](key);
	}

	std::pair<typename base_type::iterator,bool> insert(const typename base_type::value_type& val)
	{
		typename base_type::key_type & key = val.first;
    	_confirm_key(key);
		return base.insert(val);
	}

	template <class... Args>
	std::pair<typename base_type::iterator,bool> emplace (Args&&... args)
	{
		typename base_type::value_type val(args...);
		typename base_type::key_type & key = val.first;
    	_confirm_key(key);
		return base.insert(val);
	}
	template <class... Args>
	typename base_type::iterator emplace_hint (typename base_type::const_iterator position, Args&&... args)
	{
		return emplace(args...);
	}

	void erase (typename base_type::iterator position)
	{
		typename base_type::key_type key = position->first;
		base.erase(position);
		_key_map_.erase(key);
	}
	typename base_type::size_type erase (const typename base_type::key_type& k)
	{
		_update_key_map();
		auto result=base.erase(k);
		_key_map_.erase(k);
		return result;
	}
	void erase (typename base_type::iterator first, typename base_type::iterator last)
	{
		for(auto it=first;it!=last;++it)
		{
			erase(it);
		}
	}

	void swap (table_map_t<T>& other)
	{
		using std::swap;
		swap(_key_map_,other._key_map_);
		swap(_counter_,other._counter_);
		base.swap(other.base);
	}

	void clear()
	{
		_key_map_.clear();
		_counter_ = 0;
		base.clear();
	}

#endif

    // Map functions we're using as-is
#if(1)
    typename base_type::iterator begin()
    {
    	return base.begin();
    }
    typename base_type::const_iterator begin() const
    {
    	return base.begin();
    }
    typename base_type::iterator end()
    {
    	return base.end();
    }
    typename base_type::const_iterator end() const
    {
    	return base.end();
    }
    typename base_type::reverse_iterator rbegin()
    {
    	return base.rbegin();
    }
    typename base_type::const_reverse_iterator rbegin() const
    {
    	return base.rbegin();
    }
    typename base_type::reverse_iterator rend()
    {
    	return base.rend();
    }
    typename base_type::const_reverse_iterator rend() const
    {
    	return base.rend();
    }
    typename base_type::const_iterator cbegin() const noexcept
    {
    	return base.cbegin();
    }
    typename base_type::const_iterator cend() const noexcept
    {
    	return base.cend();
    }
    typename base_type::const_reverse_iterator crbegin() const noexcept
    {
    	return base.crbegin();
    }
    typename base_type::const_reverse_iterator crend() const noexcept
    {
    	return base.crend();
    }

    bool empty() const
    {
    	return base.empty();
    }
    typename base_type::size_type size() const
    {
    	return base.size();
    }
    typename base_type::size_type max_size() const
    {
    	return base.max_size();
    }
#endif
};

} // namespace brgastro

#endif // _BRG_TABLE_MAP_HPP_INCLUDED_
