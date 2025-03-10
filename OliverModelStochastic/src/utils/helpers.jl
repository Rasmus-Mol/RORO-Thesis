# Rasmus: Fraction of slot within this frame section
function calculate_slot_frame_overlap_matrix(slots, frame_positions)
    num_slots = length(slots)
    num_frames = length(frame_positions)
    overlap_matrix = zeros(Float64, num_slots, num_frames)
    
    # Assuming frames are ordered by position
    for (slot_idx, slot) in enumerate(slots)
        slot_start = slot.lcg - slot.length/2
        slot_end = slot.lcg + slot.length/2
        
        for frame_idx in 1:num_frames
            # Calculate frame boundaries
            frame_start = frame_idx == 1 ? 
                -Inf : 
                (frame_positions[frame_idx-1].position + frame_positions[frame_idx].position) / 2
            
            frame_end = frame_idx == num_frames ? 
                Inf : 
                (frame_positions[frame_idx].position + frame_positions[frame_idx+1].position) / 2
            
            # Calculate overlap
            overlap_start = max(slot_start, frame_start)
            overlap_end = min(slot_end, frame_end)
            
            if overlap_end > overlap_start
                overlap_length = overlap_end - overlap_start
                # Calculate percentage of slot that falls within this frame section
                overlap_percentage = overlap_length / slot.length
                overlap_matrix[slot_idx, frame_idx] = overlap_percentage
            end
        end
    end
    
    return overlap_matrix
end