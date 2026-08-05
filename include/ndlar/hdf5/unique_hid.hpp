#pragma once

#include <hdf5.h>

namespace ndlar::hdf5
{

/*
 * Wrapper class for managing HDF5 object identifiers (hid_t) with automatic resource management.
 */
class UniqueHID
{
public:
    using Deleter = herr_t (*)(hid_t);

    UniqueHID() :
        id_(H5I_INVALID_HID),
        closer_(nullptr)
    {
    }
    UniqueHID(hid_t id, Deleter closer) :
        id_(id),
        closer_(closer)
    {
    }

    ~UniqueHID()
    {
        reset();
    }

    UniqueHID(UniqueHID &&other) noexcept :
        id_(other.id_),
        closer_(other.closer_)
    {
        other.id_ = H5I_INVALID_HID;
    }

    UniqueHID &operator=(UniqueHID &&other) noexcept
    {
        if (this != &other)
        {
            reset();
            id_ = other.id_;
            closer_ = other.closer_;
            other.id_ = H5I_INVALID_HID;
        }
        return *this;
    }

    UniqueHID(const UniqueHID &) = delete;
    UniqueHID &operator=(const UniqueHID &) = delete;

    hid_t get() const
    {
        return id_;
    }

    operator hid_t() const
    {
        return id_;
    }

    bool is_valid() const
    {
        return id_ >= 0;
    }

    void reset(hid_t id = H5I_INVALID_HID, Deleter closer = nullptr)
    {
        if (id_ >= 0 && closer_)
        {
            closer_(id_);
        }

        id_ = id;
        closer_ = closer;
    }

private:
    hid_t id_;
    Deleter closer_;
};

} // namespace ndlar::hdf5
