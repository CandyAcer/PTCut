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

namespace MCAL {   // Mesh Cut Algorithm Libaray, 存放多面体网格切割算法的所有代码
namespace MESH {   // 存放多面体网格实现的相关代码

/*
* 属性的管理方式
*
* part 01: BasePropertyArray and PropertyArray<T>
*     我们使用动态增长的数组来存储几何元素的某一种属性, 即PropertyArray
*     几何元素的属性类型多种多样, 所以PropertyArray必须为模板类,
*     实现上, 通过对std::vector<T>适配即可达到目的.
*
* part 02: PropertyContainer
*     负责存储某种几何元素的所有属性, 例如vertex的所有属性都在同一Container中,
*     主要的作用是统一处理, 在增删几何元素时统一管理所有属性, 例如加一个元素时所有属性都push_back
*     PropertyArray<T>类型不定, 所以需要抽象出一个BasePropertyArray,
*     PropertyContainer通过若干BasePropertyArray指针管理诸多属性.
*
*     Note: PropertyContainer只负责存储, 不负责值的修改.
*
* part 03: PropertyMapBase
*     负责属性值的获取与修改, 内部有一个指向PropertyArray<T>的指针, 通过该指针修改
*     PropertyContainer修改属性值很不方便, 所以加入该类
*
*/

/*
* 几何元素的属性类型多种多样, 所以在实现上采取模板类
* 这些模板类由PropertyContainer统一管理很不方便(模板参数类型不定)
* 所以用一个基类将PropertyArray<T>抽象出来
* PropertyContainer通过基类指针即可实现统一管理
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
    // 属性名, property_name是模板类公共的数据成员, 应放在基类
    // 命名规范：'key:value', 如v:point代表vertex的坐标点属性
    std::string property_name;
};

/*
* PropertyArray<T>代表一种具体的属性
* 底层依赖vector, 通过适配完成相应功能
*
* @param T：属性值的类型
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
        // 此时属性值是默认值, 后续可通过PropertyMap修改
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
    vector_type   property_data;    // 属性对应的数组
    value_type    default_value;    // 属性的默认值, 用于初始化和reset
};


/*
* PropertyContainer负责存储并统一管理某种几何元素的若干属性.
* 一种属性对应一个vector, 多个属性对应多个vector.
* PropertyContainer主要负责对多个vector进行统一操作.
* 例如新增一个点时Container中的所有属性都需要新增一个, 重置删除同理.
*
* Note1：Container中所有vector的size和capacity是相同的.
* Note2：PropertyContainer并不对某种属性的某个值进行修改和获取, 这些功能由PropertyMap实现.
*
* @param RefClass：所在的类, 一般为PolyhedralMesh
* @param Key：几何元素(vertex, edge, face, cell)索引
*/
template <typename RefClass, typename Key>
class PropertyContainer
{
public:
    // 只需默认构造函数, 因为属性都是后续逐个添加的, 初始状态props为空
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

    // 属性数组当前的size(属性值的个数)
    std::size_t size() const { return size_; }

    // 属性数组当前的capacity
    std::size_t capacity() const { return capacity_; }

    // 属性的个数
    std::size_t n_properties() const { return props.size(); }

    // 返回所有属性名
    std::vector<std::string> properties() const
    {
        std::vector<std::string> names;
        for (std::size_t i = 0; i < props.size(); ++i)
            names.push_back(props[i]->name());
        return names;
    }

    // 删除所有属性(不只是删除props, 还要删除每个指针所指向的内存空间)
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

        // props.resize(n)之前, 需要将指针所指的内存释放, 避免内存泄露
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

    // 两个容器互换
    void swap(PropertyContainer& rhs)
    {
        this->props.swap(rhs.props);
        std::swap(this->size_, rhs.size_);
        std::swap(this->capacity_, rhs.capacity_);
    }

public:    /* 属性的添加删除获取 */

    template <typename T>
    struct GetPmapType
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type type;
    };

    // 按照名字获取第i个属性
    template <typename T>
    std::pair<typename GetPmapType<T>::type, bool>
    get(const std::string& name, std::size_t i) const
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type Pmap;

        if (props[i]->name() == name)       // 属性存在, 返回指向该属性的指针
        {
            if (PropertyArray<T>* p = dynamic_cast<PropertyArray<T>*>(props[i]))
                return std::make_pair(Pmap(p), true);
        }
        return std::make_pair(Pmap(), false); // 属性不存在, 返回nullptr
    }

    // 按照名字搜索属性
    template <typename T>
    std::pair<typename GetPmapType<T>::type, bool>
    get(const std::string& name) const
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type Pmap;
        for (std::size_t i = 0; i < props.size(); ++i)
        {
            std::pair<Pmap, bool> res = get<T>(name, i);
            if (res.second)    // 属性存在, 返回指向该属性的指针
                return res;
        }
        // 属性不存在
        return std::make_pair(Pmap(), false);
    }

    // 添加属性, name为属性名, t为属性的默认值
    template <typename T>
    std::pair<typename GetPmapType<T>::type, bool>
    add(const std::string& name, const T t = T())
    {
        typedef typename RefClass::template GetPropertyMap<Key, T>::type Pmap;
        for (std::size_t i = 0; i < props.size(); ++i)
        {
            std::pair<Pmap, bool> res = get<T>(name, i);
            if (res.second)    // 该属性已存在, 直接返回即可, false代表插入失败
            {
                res.second = false;
                return res;
            }
        }

        // 属性不存在, 新增属性, true代表插入成功
        // 注意：新增的属性所有值都是默认值, 修改由PropertyMap完成
        PropertyArray<T>* p = new PropertyArray<T>(name, t);
        p->reserve(capacity_);
        p->resize(size_);
        props.push_back(p);
        return std::make_pair(Pmap(p), true);
    }

    // 删除属性
    template <typename T>
    bool remove(typename GetPmapType<T>::type& pm)
    {
        typename std::vector<BasePropertyArray*>::iterator it = props.begin();
        for (; it != props.end(); ++it)
        {
            // 此处用了private成员, 所以PropertyContainer必须为PropertyMapBase的友元类
            if (*it == pm.parray)   // 存在该属性, 进行删除
            {
                delete* it;
                props.erase(it);
                pm.reset();
                return true;
            }
        }
        return false;              // 不存在该属性, 返回false
    }

    // 属性的值类型, 通过属性的名称查询
    // 如果属性不存在, 返回typeid(void)
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
    std::vector<BasePropertyArray*> props; // 每一个指针指向一个具体的属性
    std::size_t size_ = 0;        // size_ == props[i]->size()
    std::size_t capacity_ = 0;    // capacity_ == props[i]->capacity()
};


/*
* PropertyMapBase负责对某种属性的值进行修改或获取
* PropertyMapBase里面封装了一个指向PropertyArray<T>的指针, 通过该指针获取或修改属性值
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
    PropertyArray<T>* parray;  // 指向PropertyArray<T>的指针, 通过该指针修改属性值
};

}	// namespace MCAL::MESH
}	// namespace MCAL

#endif