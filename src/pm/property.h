// $License: GPL-3.0
// $Author: Liu Zhiqiang(from Beihang University)
// $Email: lzhq0930@gmail.com(liuzhiqiang@buaa.edu.cn)
// 
// Geometric elements can be attached many kinds of properties.
// I use dynamic array(i.e. std::vector) to store a specific property.
// In order to manage them in an uniformed and dynamic way, PropertyContainer
// and PropertyMap are introduced.
// 


#ifndef MCAL_POLYHEDRAL_MESH_PROPERTY_H
#define MCAL_POLYHEDRAL_MESH_PROPERTY_H

#include <string>
#include <vector>
#include <typeinfo>
#include <utility>
#include <assert.h>

namespace MCAL {   // Mesh Cut Algorithm Libaray, ��Ŷ����������и��㷨�����д���
namespace MESH {   // ��Ŷ���������ʵ�ֵ���ش���

/*
* ���ԵĹ���ʽ
*
* part 01: BasePropertyArray and PropertyArray<T>
*     ����ʹ�ö�̬�������������洢����Ԫ�ص�ĳһ������, ��PropertyArray
*     ����Ԫ�ص��������Ͷ��ֶ���, ����PropertyArray����Ϊģ����,
*     ʵ����, ͨ����std::vector<T>���伴�ɴﵽĿ��.
*
* part 02: PropertyContainer
*     ����洢ĳ�ּ���Ԫ�ص���������, ����vertex���������Զ���ͬһContainer��,
*     ��Ҫ��������ͳһ����, ����ɾ����Ԫ��ʱͳһ������������, �����һ��Ԫ��ʱ�������Զ�push_back
*     PropertyArray<T>���Ͳ���, ������Ҫ�����һ��BasePropertyArray,
*     PropertyContainerͨ������BasePropertyArrayָ������������.
*
*     Note: PropertyContainerֻ����洢, ������ֵ���޸�.
*
* part 03: PropertyMapBase
*     ��������ֵ�Ļ�ȡ���޸�, �ڲ���һ��ָ��PropertyArray<T>��ָ��, ͨ����ָ���޸�
*     PropertyContainer�޸�����ֵ�ܲ�����, ���Լ������
*
*/

/*
* ����Ԫ�ص��������Ͷ��ֶ���, ������ʵ���ϲ�ȡģ����
* ��Щģ������PropertyContainerͳһ����ܲ�����(ģ��������Ͳ���)
* ������һ�����ཫPropertyArray<T>�������
* PropertyContainerͨ������ָ�뼴��ʵ��ͳһ����
*/
class BasePropertyArray
{
public:
    // Default constructor.
    BasePropertyArray(const std::string& name) : property_name(name) {}

    // Destructor.
    virtual ~BasePropertyArray() {}

    // Reserve memory for n elements.
    virtual void reserve(std::size_t n) = 0;

    // Resize storage to hold n elements.
    virtual void resize(std::size_t n) = 0;

    // Free unused memory.
    virtual void shrink_to_fit() = 0;

    // Add an element in the tail of array.
    virtual void push_back() = 0;

    // Reset element to default value.
    virtual void reset(std::size_t i) = 0;

    // Let two elements swap their storage place.
    virtual void swap(std::size_t i0, std::size_t i1) = 0;

    // Return a deep copy of self.
    virtual BasePropertyArray* clone() const = 0;

    // Return a empty copy of self.
    virtual BasePropertyArray* empty_clone() const = 0;

    // Return the type_info of the property.
    virtual const std::type_info& type() const = 0;

    // Return the name of the property.
    const std::string& name() const { return property_name; }

    // Judge whether two property are equal.
    bool is_same(const BasePropertyArray& rhs)
    {
        return (name() == rhs.name() && type() == rhs.type());
    }

protected:
    // ������, property_name��ģ���๫�������ݳ�Ա, Ӧ���ڻ���
    // �����淶��'key:value', ��v:point����vertex�����������
    std::string property_name;
};

/*
* PropertyArray<T>����һ�־��������
* �ײ�����vector, ͨ�����������Ӧ����
*
* @param T������ֵ������
*/
template <typename T>
class PropertyArray : public BasePropertyArray
{
public:
    typedef T                                       value_type;
    typedef std::vector<value_type>                 vector_type;
    typedef typename vector_type::reference         reference;
    typedef typename vector_type::const_reference   const_reference;
    typedef typename vector_type::iterator          iterator;
    typedef typename vector_type::const_iterator    const_iterator;

    PropertyArray(const std::string& name, T t = T())
        : BasePropertyArray(name), default_value(t)
    {}

public:

    // Reserve memory for n elements.
    virtual void reserve(std::size_t n)
    {
        property_data.reserve(n);
    }

    // Resize storage to hold n elements.
    virtual void resize(std::size_t n)
    {
        property_data.resize(n);
    }

    // Free unused memory.
    virtual void shrink_to_fit()
    {
        vector_type(property_data).swap(property_data);
    }

    // Add an element in the tail of array.
    virtual void push_back()
    {
        // ��ʱ����ֵ��Ĭ��ֵ, ������ͨ��PropertyMap�޸�
        property_data.push_back(default_value);
    }

    // Reset element to default value.
    virtual void reset(std::size_t i)
    {
        property_data[i] = default_value;
    }

    // Let two elements swap their storage place.
    virtual void swap(std::size_t i0, std::size_t i1)
    {
        T d(property_data[i0]);
        property_data[i0] = property_data[i1];
        property_data[i1] = d;
    }

    // Return a deep copy of self.
    virtual BasePropertyArray* clone() const
    {
        PropertyArray<T>* p = 
            new PropertyArray<T>(this->property_name, this->default_value);
        p->property_data = property_data;
        return p;
    }

    // Return a empty copy of self.
    virtual BasePropertyArray* empty_clone() const
    {
        PropertyArray<T>* p = 
            new PropertyArray<T>(this->property_name, this->default_value);
        return p;
    }

    // Return the type_info of the property.
    virtual const std::type_info& type() const { return typeid(T); }

public:
    // Get pointer to array (does not work for T==bool)
    const T* data() const
    {
        return &property_data[0];
    }

    // Access the i'th element.
    reference operator[](std::size_t i)
    {
        assert(i < property_data.size());
        return property_data[i];
    }

    // Const access to the i'th element.
    const_reference operator[](std::size_t i) const
    {
        assert(i < property_data.size());
        return property_data[i];
    }

    iterator begin() { return property_data.begin(); }
    iterator end() { return property_data.end(); }
    const_iterator begin() const { return property_data.begin(); }
    const_iterator end() const { return property_data.end(); }

private:
    vector_type   property_data;    // ���Զ�Ӧ������
    value_type    default_value;    // ���Ե�Ĭ��ֵ, ���ڳ�ʼ����reset
};


/*
* PropertyContainer����洢��ͳһ����ĳ�ּ���Ԫ�ص���������.
* һ�����Զ�Ӧһ��vector, ������Զ�Ӧ���vector.
* PropertyContainer��Ҫ����Զ��vector����ͳһ����.
* ��������һ����ʱContainer�е��������Զ���Ҫ����һ��, ����ɾ��ͬ��.
*
* Note1��Container������vector��size��capacity����ͬ��.
* Note2��PropertyContainer������ĳ�����Ե�ĳ��ֵ�����޸ĺͻ�ȡ, ��Щ������PropertyMapʵ��.
*
* @param RefClass�����ڵ���, һ��ΪPolyhedralMesh
* @param Key������Ԫ��(vertex, edge, face, cell)����
*/
template <typename RefClass, typename Key>
class PropertyContainer
{
public:
    // ֻ��Ĭ�Ϲ��캯��, ��Ϊ���Զ��Ǻ��������ӵ�, ��ʼ״̬propsΪ��
    PropertyContainer() = default;

    virtual ~PropertyContainer() { clear(); }

    PropertyContainer(const PropertyContainer& rhs) { operator=(rhs); }

    PropertyContainer& operator=(const PropertyContainer& rhs)
    {
        if (this != &rhs)
        {
            clear();
            props.resize(rhs.n_properties());
            size_ = rhs.size();
            capacity_ = rhs.capacity();
            for (std::size_t i = 0; i < props.size(); ++i)
                props[i] = rhs.props[i]->clone();
        }
        return *this;
    }

    // �������鵱ǰ��size(����ֵ�ĸ���)
    std::size_t size() const { return size_; }

    // �������鵱ǰ��capacity
    std::size_t capacity() const { return capacity_; }

    // ���Եĸ���
    std::size_t n_properties() const { return props.size(); }

    // ��������������
    std::vector<std::string> properties() const
    {
        std::vector<std::string> names;
        for (std::size_t i = 0; i < props.size(); ++i)
            names.push_back(props[i]->name());
        return names;
    }

    // ɾ����������(��ֻ��ɾ��props, ��Ҫɾ��ÿ��ָ����ָ����ڴ�ռ�)
    void clear()
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            delete props[i];
        props.clear();
        size_ = 0;
    }

    // Reserve memory for n entries in all arrays.
    void reserve(std::size_t n)
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            props[i]->reserve(n);
        capacity_ = std::max(n, capacity_);
    }

    // Resize all arrays to size n.
    void resize(std::size_t n)
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            props[i]->resize(n);
        size_ = n;
    }

    // Resize the number of properties to n, deleting all other properties.
    void resize_property_array(std::size_t n)
    {
        if (props.size() <= n)
            return;

        // props.resize(n)֮ǰ, ��Ҫ��ָ����ָ���ڴ��ͷ�, �����ڴ�й¶
        for (std::size_t i = n; i < props.size(); ++i)
            delete props[i];
        props.resize(n);
    }

    // Free unused space in all arrays.
    void shrink_to_fit()
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            props[i]->shrink_to_fit();
        capacity_ = size_;
    }

    // Add a new element to each vector.
    void push_back()
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            props[i]->push_back();
        ++size_;
        capacity_ = std::max(size_, capacity_);
    }

    // Reset element to its default property values.
    void reset(std::size_t idx)
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            props[i]->reset(idx);
    }

    // Swap elements i0 and i1 in all arrays.
    void swap(std::size_t i0, std::size_t i1)
    {
        for (std::size_t i = 0; i < props.size(); ++i)
            props[i]->swap(i0, i1);
    }

    // ������������
    void swap(PropertyContainer& rhs)
    {
        this->props.swap(rhs.props);
        std::swap(this->size_, rhs.size_);
        std::swap(this->capacity_, rhs.capacity_);
    }

public:    /* ���Ե����ɾ����ȡ */

    template <typename T>
    struct GetPmapType
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type type;
    };

    // �������ֻ�ȡ��i������
    template <typename T>
    std::pair<typename GetPmapType<T>::type, bool>
    get(const std::string& name, std::size_t i) const
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type Pmap;

        if (props[i]->name() == name)       // ���Դ���, ����ָ������Ե�ָ��
        {
            if (PropertyArray<T>* p = dynamic_cast<PropertyArray<T>*>(props[i]))
                return std::make_pair(Pmap(p), true);
        }
        return std::make_pair(Pmap(), false); // ���Բ�����, ����nullptr
    }

    // ����������������
    template <typename T>
    std::pair<typename GetPmapType<T>::type, bool>
    get(const std::string& name) const
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type Pmap;
        for (std::size_t i = 0; i < props.size(); ++i)
        {
            std::pair<Pmap, bool> res = get<T>(name, i);
            if (res.second)    // ���Դ���, ����ָ������Ե�ָ��
                return res;
        }
        // ���Բ�����
        return std::make_pair(Pmap(), false);
    }

    // �������, nameΪ������, tΪ���Ե�Ĭ��ֵ
    template <typename T>
    std::pair<typename GetPmapType<T>::type, bool>
    add(const std::string& name, const T t = T())
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type Pmap;
        for (std::size_t i = 0; i < props.size(); ++i)
        {
            std::pair<Pmap, bool> res = get<T>(name, i);
            if (res.second)    // �������Ѵ���, ֱ�ӷ��ؼ���, false�������ʧ��
            {
                res.second = false;
                return res;
            }
        }

        // ���Բ�����, ��������, true�������ɹ�
        // ע�⣺��������������ֵ����Ĭ��ֵ, �޸���PropertyMap���
        PropertyArray<T>* p = new PropertyArray<T>(name, t);
        p->reserve(capacity_);
        p->resize(size_);
        props.push_back(p);
        return std::make_pair(Pmap(p), true);
    }

    // ɾ������
    template <typename T>
    bool remove(typename GetPmapType<T>::type& pm)
    {
        typename std::vector<BasePropertyArray*>::iterator it = props.begin();
        for (; it != props.end(); ++it)
        {
            // �˴�����private��Ա, ����PropertyContainer����ΪPropertyMapBase����Ԫ��
            if (*it == pm.parray)   // ���ڸ�����, ����ɾ��
            {
                delete* it;
                props.erase(it);
                pm.reset();
                return true;
            }
        }
        return false;              // �����ڸ�����, ����false
    }

    // ���Ե�ֵ����, ͨ�����Ե����Ʋ�ѯ
    // ������Բ�����, ����typeid(void)
    const std::type_info&
    get_type(const std::string& name) const
    {
        for (std::size_t i = 0; i < props.size(); ++i)
        {
            if (props[i]->name() == name)
                return props[i]->type();
        }
        return typeid(void);
    }

private:
    std::vector<BasePropertyArray*> props; // ÿһ��ָ��ָ��һ�����������
    std::size_t size_ = 0;        // size_ == props[i]->size()
    std::size_t capacity_ = 0;    // capacity_ == props[i]->capacity()
};


/*
* PropertyMapBase�����ĳ�����Ե�ֵ�����޸Ļ��ȡ
* PropertyMapBase�����װ��һ��ָ��PropertyArray<T>��ָ��, ͨ����ָ���ȡ���޸�����ֵ
*
* @param I: The key type of the property map.
* @param T: The value type of the property map.
*
*/
template <typename I, typename T>
class PropertyMapBase
{
public:
    template <typename RefClass, typename Key>
    friend class PropertyContainer;

    typedef I key_type;
    typedef T value_type;

    typedef typename PropertyArray<value_type>::reference          reference;
    typedef typename PropertyArray<value_type>::const_reference    const_reference;
    typedef typename PropertyArray<value_type>::iterator           iterator;
    typedef typename PropertyArray<value_type>::const_iterator     const_iterator;

public:
    PropertyMapBase(PropertyArray<T>* p = nullptr) : parray(p) {}

    PropertyMapBase(const PropertyMapBase& pm) : parray(pm.parray) {}

    PropertyMapBase& operator=(const PropertyMapBase& pm)
    {
        parray = pm.parray;
        return *this;
    }

    void reset()
    {
        parray = nullptr;
    }

    const T* data() const
    {
        assert(parray != nullptr);
        return parray->data();
    }

    reference operator[](const key_type& k)
    {
        assert(parray != nullptr);
        return (*parray)[k];
    }

    const_reference operator[](const key_type& k) const
    {
        assert(parray != nullptr);
        return (*parray)[k];
    }

    iterator begin() { return parray->begin(); }
    iterator end() { return parray->end(); }
    const_iterator begin() const { return parray->begin(); }
    const_iterator end() const { return parray->end(); }

private:
    PropertyArray<T>* parray;  // ָ��PropertyArray<T>��ָ��, ͨ����ָ���޸�����ֵ
};

}	// namespace MCAL::MESH
}	// namespace MCAL

#endif