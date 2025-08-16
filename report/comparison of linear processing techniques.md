**Bảng 1.1: Bảng so sánh hiệu suất các kỹ thuật kết hợp thu trong UL Massive MIMO**

| Tiêu chí             | LMMSE                          | ZF                              | MRC                              |
|-----------------------|--------------------------------|----------------------------------|----------------------------------|
| Tăng số lượng anten   | Hiệu suất tăng mạnh            | Tăng nhưng cần (M >> K)          | Có cải thiện nhưng giới hạn      |
| SNR thấp              | Tốt, cân bằng nhiễu & can nhiễu| Kém, khuếch đại nhiễu           | Ổn định nhưng còn can nhiễu      |
| SNR cao               | Tiệm cận ZF                   | Rất tốt, triệt tiêu can nhiễu    | Hiệu suất giảm do can nhiễu      |
| Khử can nhiễu         | Có một phần                   | Triệt tiêu hoàn toàn             | Không khử được                   |
| Độ phức tạp           | Cao                           | Trung bình                       | Thấp                             |

**Bảng 1.2: Bảng so sánh hiệu suất các kỹ thuật tiền mã hóa trong DL Massive MIMO**

| Tiêu chí             | LMMSE                | ZF                          | RZF                                | MRT                                |
|-----------------------|----------------------|-----------------------------|-------------------------------------|-------------------------------------|
| Tăng số lượng anten   | Hiệu suất tăng mạnh  | Tăng hiệu suất khi M >> K   | Tối ưu so với ZF                   | Có cải thiện nhưng còn giới hạn     |
| SNR thấp              | Rất tốt              | Kém                         | Tốt                                 | Tốt                                 |
| SNR cao               | Gần giống ZF         | Rất tốt                     | Tiệm cận ZF                         | Can nhiễu lớn, giảm hiệu suất       |
| Khử can nhiễu         | Có                   | Có                          | Có (phụ thuộc hệ số điều chuẩn)    | Không                               |
| Độ phức tạp           | Cao                  | Trung bình                  | Trung bình - cao                    | Thấp                                |

**Bảng 1.3: Bảng so sánh hiệu suất kỹ thuật tiền mã hóa ZF trong DL Massive MIMO**

| Tiêu chí             | ZF - Phân bổ công suất đều                 | ZF - Phân bổ công suất tối ưu            |
|-----------------------|---------------------------|-------------------------|
| Tăng số lượng anten   | Kém                       | Tốt                     |
| SNR thấp              | Kém                       | Tốt                     |
| SNR cao               | Tương tự như phân bổ công suất đều                  | Tốt                     |
| Hiệu suất tổng thể    | Giới hạn ở SNR thấp       | Tối ưu                  |
