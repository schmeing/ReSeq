#ifndef BOUNDEDWORKQUEUE_H
#define BOUNDEDWORKQUEUE_H

#include <array>
#include <atomic>
#include <cstddef>
#include <mutex>
#include <vector>

namespace reseq {

struct QueueSlot {
    std::atomic<uint32_t> current_param{0};
    std::atomic<uint32_t> finished_count{0};
    uint32_t total_params{0};
    std::atomic<bool> published{false};
};

template <size_t MaxSlots> class BoundedWorkQueue {
    std::array<QueueSlot, MaxSlots> slots_;
    std::vector<size_t> free_indices_;
    mutable std::mutex mutex_;

  public:
    BoundedWorkQueue() {
        free_indices_.reserve(MaxSlots);
        for (size_t i = MaxSlots; i--;) {
            free_indices_.push_back(i);
        }
    }

    size_t try_acquire() {
        std::lock_guard lock(mutex_);
        if (free_indices_.empty()) {
            return SIZE_MAX;
        }
        size_t idx = free_indices_.back();
        free_indices_.pop_back();
        return idx;
    }

    void publish(size_t idx, uint32_t total) {
        slots_[idx].finished_count.store(0, std::memory_order_relaxed);
        slots_[idx].total_params = total;
        slots_[idx].current_param.store(0, std::memory_order_relaxed);
        // published must be set LAST with release semantics so consumers
        // reading published.load(acquire) see all preceding writes.
        slots_[idx].published.store(true, std::memory_order_release);
    }

    void release_empty(size_t idx) {
        std::lock_guard lock(mutex_);
        slots_[idx].published = false;
        free_indices_.push_back(idx);
    }

    void release(size_t idx) {
        std::lock_guard lock(mutex_);
        slots_[idx].published = false;
        free_indices_.push_back(idx);
    }

    QueueSlot& slot(size_t idx) { return slots_[idx]; }
    const QueueSlot& slot(size_t idx) const { return slots_[idx]; }
    static constexpr size_t max_slots() { return MaxSlots; }
};

} // namespace reseq

#endif // BOUNDEDWORKQUEUE_H
