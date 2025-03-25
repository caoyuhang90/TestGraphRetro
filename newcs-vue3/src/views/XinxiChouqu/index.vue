<template>
  <el-container style="height: 100vh;">
    <!-- 顶部导航栏 -->
    <el-header class="header">
      化学反应知识抽取
    </el-header>

    <!-- 主体内容 -->
    <el-main style="padding: 20px; display: flex; justify-content: space-between; flex-wrap: wrap;">
      <el-card class="module" style="width: 100%">
        <div>
          <el-text class="mx-1" size="large">请选择需要抽取的内容（文本、图片、文献）</el-text>
        </div>
        <el-radio-group v-model="inputMethod">
          <el-radio label="string">文本抽取</el-radio>
          <el-radio label="image">分子结构识别</el-radio>
          <el-radio label="reaction">反应抽取</el-radio>
          <el-radio label="article">文献抽取</el-radio>
        </el-radio-group>

        <!-- 如果选择了输入字符串方式 -->
        <div v-if="inputMethod === 'string'">
          <div>
            <el-form label-position="top" style="margin-top: 20px;">
              <el-form-item label="请输入一段描述化学反应的文本信息：">
                <!-- 输入1 -->
                <el-input v-model="inputText" autosize type="textarea"
                  placeholder="例如：Finally, porphyrin 4 was quaternized by reaction with a large excessof methyl triﬂate in hot DMF to yield the desired Fe-o-TMA as atriﬂate salt."></el-input>
              </el-form-item>
              <el-form-item>
                <el-button type="primary" @click="submitTextExtraction">开始抽取信息</el-button>
              </el-form-item>
            </el-form>
          </div>

          <div>
            <el-text class="mx-1" size="large">示例数据</el-text>
            <el-table :data="exampleData1" style="width: 100%" height="auto">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column prop="reactants" label="底物Reactants" width="200" align="center" />
              <el-table-column prop="product" label="产物Product" width="200" align="center" />
              <el-table-column prop="solvent" label="溶剂Solvent" width="200" align="center" />
              <el-table-column prop="temperature" label="温度Temperature" width="200" align="center" />
              <el-table-column prop="catalyst" label="催化剂Catalyst" width="200" align="center" />
              <el-table-column prop="time" label="时间Time" width="200" align="center" />
              <el-table-column prop="yield" label="产量Yield" width="200" align="center" />
            </el-table>
          </div>

          <div v-if="textExtractionResult">
            <el-text class="mx-1" size="large">预测结果</el-text>
            <el-table :data="textExtractionResult.data" style="width: 100%" height="auto">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column prop="reactants" label="底物Reactants" width="200" align="center" />
              <el-table-column prop="product" label="产物Product" width="200" align="center" />
              <el-table-column prop="solvent" label="溶剂Solvent" width="200" align="center" />
              <el-table-column prop="temperature" label="温度Temperature" width="200" align="center" />
              <el-table-column prop="catalyst" label="催化剂Catalyst" width="200" align="center" />
              <el-table-column prop="time" label="时间Time" width="200" align="center" />
              <el-table-column prop="yield" label="产量Yield" width="200" align="center" />
            </el-table>

          </div>
        </div>


        <div v-if="inputMethod === 'image'">
          <div>
            <el-form label-position="top" style="margin-top: 20px;">
              <el-form-item label="请上传一张分子结构图(PNG、JPG格式)：">
                <el-upload class="upload-demo" drag @change="handleImageChange" accept=".jpg, .png" :auto-upload="false"
                  multiple list-type="picture-card">
                  <i class="el-icon-upload"></i>
                  <div class="el-upload__text">点击或拖拽上传</div>
                </el-upload>
              </el-form-item>
              <el-form-item>
                <el-button type="primary" @click="submitSmilesExtraction" :loading="uploading">开始抽取信息</el-button>
              </el-form-item>
            </el-form>
          </div>
          <div>
            <el-text class="mx-1" size="large">示例数据</el-text>
            <el-table :data="exampleData2" style="width: 100%" height="auto">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column label="分子结构图" width="600" align="center">
                <template #default="scope">
                  <el-image style="width: 100px; height: 100px" :src="scope.row.image"></el-image>
                </template>
              </el-table-column>
              <el-table-column prop="smiles" label="SMILES表达式" width="600" align="center" />
            </el-table>
          </div>
          <div v-if=smilesExtractionResult>
            <el-text class="mx-1" size="large">预测结果</el-text>
            <el-table :data="smilesExtractionResult.data" style="width: 100%" height="auto">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column label="分子结构图" width="600" align="center">
                <template #default="scope">
                  <el-image style="width: 100px; height: 100px" :src="scope.row.image"></el-image>
                </template>
              </el-table-column>
              <el-table-column prop="smiles" label="SMILES表达式" width="600" align="center" />
            </el-table>
          </div>
        </div>
        <!-- 如果选择了上传文件方式 -->

        <div v-if="inputMethod === 'reaction'">
          <div>
            <el-form label-position="top" style="margin-top: 20px;">
              <el-form-item label="请上传一张化学反应方程的图片(PNG、JPG格式)：">
                <el-upload class="upload-demo" drag @change="handleImageChange" accept=".jpg, .png" :auto-upload="false"
                  multiple list-type="picture-card">
                  <i class="el-icon-upload"></i>
                  <div class="el-upload__text">点击或拖拽上传</div>
                </el-upload>
              </el-form-item>
              <el-form-item>
                <el-button type="primary" @click="submitReactionExtraction">开始抽取信息</el-button>
              </el-form-item>
            </el-form>
          </div>
          <div>
            <el-text class="mx-1" size="large">示例数据</el-text>
            <el-table :data="exampleData3" style="width: 100%" height="300">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column label="底物分子结构图" width="400" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="reactants" label="反应底物" width="250" align="center" />
              <el-table-column prop="conditions" label="反应条件" width="250" align="center" />
              <el-table-column label="产物分子结构图" width="400" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="products" label="产物" width="250" align="center" />
            </el-table>
          </div>
          <div v-if="reactionExtractionResult">
            <el-text class="mx-1" size="large">预测结果</el-text>
            <el-table :data="reactionExtractionResult.data" style="width: 100%" height="300">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column label="底物分子结构图" width="400" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="reactants" label="反应底物" width="250" align="center" />
              <el-table-column prop="conditions" label="反应条件" width="250" align="center" />
              <el-table-column label="产物分子结构图" width="400" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="products" label="产物" width="250" align="center" />
            </el-table>
          </div>

        </div>


        <div v-if="inputMethod === 'article'">
          <!-- 如果选择了上传文献方式 -->
          <div>
            <el-form label-position="top" style="margin-top: 20px;">
              <el-form-item label="请上传与化学相关的文献（PDF格式）:">
                <el-upload class="upload-demo" drag :before-upload="() => false"
                  @change="handleArticleChange($event, 'smilesFile')" accept=".pdf">
                  <i class="el-icon-upload"></i>
                  <div class="el-upload__text">点击或拖拽文件上传</div>
                </el-upload>
              </el-form-item>
              <el-form-item>
                <el-button type="primary" @click="submitArticleExtraction">开始抽取信息</el-button>
              </el-form-item>
            </el-form>
          </div>

          <div>
            <el-text class="mx-1" size="large">示例数据</el-text>
            <el-table :data="exampleData4" style="width: 100%" height="auto">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column label="底物分子结构图" width="200" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="reactants" label="底物" width="200" align="center" />
              <el-table-column label="产物分子结构图" width="200" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="product" label="产物" width="200" align="center" />
              <el-table-column prop="solvent" label="溶剂" width="200" align="center" />
              <el-table-column prop="temperature" label="温度" width="200" align="center" />
              <el-table-column prop="catalyst" label="催化剂" width="200" align="center" />
              <el-table-column prop="time" label="时间" width="200" align="center" />
              <el-table-column prop="yield" label="产量" width="200" align="center" />
            </el-table>
          </div>

          <div v-if="articleExtractionResult">
            <el-text class="mx-1" size="large">示例数据</el-text>
            <el-table :data="articleExtractionResult.data" style="width: 100%" height="auto">
              <el-table-column fixed prop="id" label="序号" width="100" align="center" />
              <el-table-column label="底物分子结构图" width="200" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="reactants" label="底物" width="200" align="center" />
              <el-table-column label="产物分子结构图" width="200" align="center">
                <template #default="scope">
                  <el-image v-for="(item, index) in scope.row.productImages" :key="index" :src="item"
                    style="width: 100px; height: 100px" :preview-src-list="[item]">
                  </el-image>
                </template>
              </el-table-column>
              <el-table-column prop="product" label="产物" width="200" align="center" />
              <el-table-column prop="solvent" label="溶剂" width="200" align="center" />
              <el-table-column prop="temperature" label="温度" width="200" align="center" />
              <el-table-column prop="catalyst" label="催化剂" width="200" align="center" />
              <el-table-column prop="time" label="时间" width="200" align="center" />
              <el-table-column prop="yield" label="产量" width="200" align="center" />
            </el-table>

          </div>
        </div>


      </el-card>


    </el-main>
  </el-container>
</template>


<script>
import { ElMessage } from 'element-plus';
import axios from 'axios';

let tableData;
export default {


  data() {
    return {
      inputText: '',
      activeTab: 'createTask',
      inputMethod: 'string', // 默认为输入字符串的方式


      // 识别文本
      textExtractionResult: null,
      // 识别SMILES
      smilesExtractionResult: null,
      // 识别反应方程
      reactionExtractionResult: null,
      // 识别文献
      articleExtractionResult: null,

      selectedFile: null,  // 保存选中的文件
      uploading: false,

      exampleData1: [
        {
          id: '1',
          reactants: 'reactants',
          product: 'product',
          solvent: 'solvent',
          temperature: 'temperature',
          catalyst: 'catalyst',
          time: 'time',
          yield: 'yield',
        },
        {
          id: '2',
          reactants: 'reactants',
          product: 'product',
          solvent: 'solvent',
          temperature: 'temperature',
          catalyst: 'catalyst',
          time: 'time',
          yield: 'yield',
        }
      ],
      exampleData2: [
        {
          id: '1',
          image: 'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg',
          smiles: '示例SMILES'
        }, {
          id: '2',
          image: 'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg',
          smiles: '示例SMILES'
        }

      ],
      exampleData3: [
        {
          id: '1',
          reactants: 'reactants',
          conditions: 'conditions',
          products: 'products',
          reactantImages: [
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100',
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg'],
          productImages: [
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100',
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg']
        }, {
          id: '2',
          reactants: 'reactants',
          conditions: 'conditions',
          products: 'products',
          reactantImages: [
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100',
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg'],
          productImages: [
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100',
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg']
        }

      ],
      exampleData4: [
        {
          id: '1',
          reactants: 'reactants',
          product: 'product',
          solvent: 'solvent',
          temperature: 'temperature',
          catalyst: 'catalyst',
          time: 'time',
          yield: 'yield',
          reactantImages: [
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100',
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg'],
          productImages: [
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100',
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg']
        },
        {
          id: '2',
          reactants: 'reactants',
          product: 'product',
          solvent: 'solvent',
          temperature: 'temperature',
          catalyst: 'catalyst',
          time: 'time',
          yield: 'yield',
          reactantImages: [
            'https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg',
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100'
          ],
          productImages: ['https://fuss10.elemecdn.com/e/5d/4a731a90594a4af544c0c25941171jpeg.jpeg',
            'https://fuss10.elemecdn.com/3/63/4e7f3a15429bfda99bce42a18cdd1jpeg.jpeg?imageMogr2/thumbnail/360x360/format/webp/quality/100'
          ]
        }
      ]
    };
  },
  methods: {
    async submitTextExtraction() {
      try {
        // 构建请求体
        const requestBody = {
          inputText: this.inputText
        };

        // 发送POST请求到后端
        const response = await axios.post('http://localhost:5006/api/extraction/text', requestBody);

        // 打印响应，确认数据结构
        console.log(response.data);
        this.textExtractionResult = response.data;

        // 解析返回的结果并显示在前端
      } catch (error) {
        console.error("Error submitting form:", error);
      }
    },

    // 获取用户上传的文件
    handleImageChange(file) {
      // file 是用户上传的文件对象
      this.selectedFile = file.raw;  // 通过 file.raw 获取原始文件对象
    },

    async submitSmilesExtraction() {
      if (!this.selectedFile) {
        console.error("No file selected");
        return;
      }

      const formData = new FormData();
      formData.append("file", this.selectedFile);  // 使用 'file' 作为 key，FastAPI 需要和后端一致

      try {
        this.uploading = true;
        const response = await axios.post('http://localhost:5006/api/extraction/smiles', formData, {
          headers: { 'Content-Type': 'multipart/form-data' }
        });
        console.log(response.data);
        this.smilesExtractionResult = response.data;
      } catch (error) {
        console.error("Error submitting form:", error);
      } finally {
        this.uploading = false;
      }
    },

    async submitReactionExtraction() {
      if (!this.selectedFile) {
        console.error("No file selected");
        return;
      }

      const formData = new FormData();
      formData.append("file", this.selectedFile);  // 将文件添加到 FormData

      try {
        const response = await axios.post('http://localhost:5006/api/extraction/reaction', formData, {
          headers: { 'Content-Type': 'multipart/form-data' }
        });
        console.log(response.data);
        this.reactionExtractionResult = response.data;
      } catch (error) {
        console.error("Error submitting form:", error);
      }
    },

    async submitArticleExtraction() {
      const formData = new FormData();
      formData.append("file", this.selectedFile);  // 上传的PDF文件

      const response = await axios.post('http://localhost:5006/api/extraction/article', formData, {
        headers: { 'Content-Type': 'multipart/form-data' }
      });

      console.log(response.data);
      this.articleExtractionResult = response.data;
    }
  },
};

</script>


<style scoped>
.header {
  background-color: rgb(1,70,138);
  color: white;
  text-align: center;
  padding: 0px 0;
}

.el-tabs {
  background-color: #ffffff;
  margin-bottom: 20px;
}

.module {
  padding: 20px;
  box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
  border-radius: 8px;
}

.el-upload__text {
  font-size: 12px;
  text-align: center;
}

.el-button {
  width: 100%;
}
</style>
